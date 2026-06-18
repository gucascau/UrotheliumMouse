#!/usr/bin/env python3
"""
Script 05: scVI integration of bladder urothelium per-sample QC h5ad files

Workflow
--------
1. Load all *.h5ad from qc_h5ad/ (produced by 03_export_h5ad_bladder.R).
   Samples without a raw counts layer are skipped automatically.
2. Concatenate into a single AnnData (outer join on genes).
3. Normalize + log1p on X for HVG selection; keep raw counts in layers.
4. Select 3000 HVGs (batch-aware, flavour="seurat_v3").
5. Subset to HVGs and train scVI (batch_key = gsm_id,
   categorical covariate = technology + sample_id).
6. Extract latent Z, build neighbours, run UMAP + Leiden clustering.
7. Optionally fine-tune scANVI if cell_type_original labels are available.
8. Save integrated AnnData + model directories.

Outputs (under integration_output/)
-------------------------------------
  bladder_scvi_integrated.h5ad  — integrated AnnData
  bladder_scvi_model/           — trained scVI model
  bladder_scanvi_model/         — trained scANVI model (if run)
  bladder_scvi_plots/           — UMAP PDFs + ELBO curve

Requirements:
  scvi-tools>=1.1  scanpy>=1.9  anndata  torch  leidenalg  igraph  matplotlib
"""

import os
import glob
import logging
import numpy as np
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
import scvi
import torch
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR     = "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/BladderUrothelium"
QC_H5AD_DIR  = os.path.join(DATA_DIR, "qc_h5ad")
OUT_DIR      = os.path.join(DATA_DIR, "integration_output")
PLOT_DIR     = os.path.join(OUT_DIR, "bladder_scvi_plots")
SCVI_MODEL   = os.path.join(OUT_DIR, "bladder_scvi_model")
SCANVI_MODEL = os.path.join(OUT_DIR, "bladder_scanvi_model")
H5AD_OUT     = os.path.join(OUT_DIR, "bladder_scvi_integrated.h5ad")

for d in (OUT_DIR, PLOT_DIR):
    os.makedirs(d, exist_ok=True)

# ── Parameters ─────────────────────────────────────────────────────────────────
N_HVG        = 3000
N_LATENT     = 20
N_LAYERS     = 2
N_HIDDEN     = 256
N_EPOCHS     = 400
BATCH_SIZE   = 256
LR           = 1e-3
LEIDEN_RES   = 0.5
N_NEIGHBOURS = 15

BATCH_KEY    = "gsm_id"
CAT_COVS     = ["technology", "sample_id"]

CELL_TYPE_KEY         = "cell_type_original"
UNLABELED_CAT         = "Unknown"
SCANVI_LABEL_MIN_FRAC = 0.01

EXCLUDE_SAMPLES = set()   # all bladder samples are included

# ── Logging ────────────────────────────────────────────────────────────────────
logging.basicConfig(
    level   = logging.INFO,
    format  = "%(asctime)s  %(levelname)s  %(message)s",
    datefmt = "%H:%M:%S",
)
log = logging.getLogger(__name__)
scvi.settings.verbosity = 2

# ── Device ─────────────────────────────────────────────────────────────────────
device = "cuda" if torch.cuda.is_available() else "cpu"
log.info("PyTorch device: %s", device)
if device == "cuda":
    log.info("  GPU: %s", torch.cuda.get_device_name(0))
scvi.settings.dl_num_workers = int(os.environ.get("SLURM_CPUS_PER_TASK", 4))


###############################################################################
# STEP 1: Load per-sample h5ad files
###############################################################################

h5ad_files = sorted(glob.glob(os.path.join(QC_H5AD_DIR, "*.h5ad")))
log.info("Found %d h5ad files in %s", len(h5ad_files), QC_H5AD_DIR)

adatas  = []
skipped = []

for path in h5ad_files:
    sid = os.path.splitext(os.path.basename(path))[0]

    if sid in EXCLUDE_SAMPLES:
        log.info("  [skip – excluded] %s", sid)
        skipped.append((sid, "excluded"))
        continue

    a = sc.read_h5ad(path)

    # Promote integer X to counts layer if not already split
    if "counts" not in a.layers:
        _x    = a.X
        _data = _x.data if sp.issparse(_x) else _x.ravel()
        if len(_data) > 0 and np.all(_data == np.floor(_data)):
            a.layers["counts"] = a.X.copy()
            log.info("  [counts from X] %s", sid)
        else:
            log.info("  [skip – no counts] %s  (%d cells)", sid, a.n_obs)
            skipped.append((sid, "no counts layer"))
            continue

    # Fall back to sample_id if batch key is missing
    if BATCH_KEY not in a.obs.columns or a.obs[BATCH_KEY].isna().all():
        log.warning("  '%s' missing for %s — using sample_id as batch",
                    BATCH_KEY, sid)
        a.obs[BATCH_KEY] = sid

    log.info("  [loaded] %-28s  %6d cells × %6d genes", sid, a.n_obs, a.n_vars)
    adatas.append(a)

log.info("\nSamples loaded : %d", len(adatas))
log.info("Samples skipped: %d", len(skipped))
for s, reason in skipped:
    log.info("  %-30s  (%s)", s, reason)

if len(adatas) == 0:
    raise RuntimeError("No samples with raw counts found in %s" % QC_H5AD_DIR)


###############################################################################
# STEP 2: Concatenate
###############################################################################

log.info("\nConcatenating ...")
adata = ad.concat(
    adatas,
    join         = "outer",
    merge        = "same",
    index_unique = "-",
)

# Fill NaN introduced by outer join with 0
for layer in list(adata.layers.keys()):
    if sp.issparse(adata.layers[layer]):
        adata.layers[layer].data = np.nan_to_num(adata.layers[layer].data)
    else:
        adata.layers[layer] = np.nan_to_num(adata.layers[layer])

del adatas
log.info("Concatenated: %d cells × %d genes", adata.n_obs, adata.n_vars)
log.info("Samples: %s", sorted(adata.obs[BATCH_KEY].unique()))


###############################################################################
# STEP 3: Normalize + select HVGs (batch-aware)
###############################################################################

adata.layers["counts"] = adata.layers["counts"].copy()

log.info("Normalizing for HVG selection ...")
adata.X = adata.layers["counts"].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers["log_norm"] = adata.X.copy()

log.info("Selecting %d HVGs (batch_key=%s) ...", N_HVG, BATCH_KEY)
sc.pp.highly_variable_genes(
    adata,
    n_top_genes = N_HVG,
    batch_key   = BATCH_KEY,
    flavor      = "seurat_v3",
    layer       = "counts",
    subset      = False,
)
log.info("  HVGs selected: %d", adata.var["highly_variable"].sum())


###############################################################################
# STEP 4: Subset to HVGs + set up scVI
###############################################################################

adata_hvg = adata[:, adata.var["highly_variable"]].copy()
log.info("After HVG subset: %d cells × %d genes", adata_hvg.n_obs, adata_hvg.n_vars)

adata_hvg.X = adata_hvg.layers["counts"]

log.info("Setting up scVI ...")
cat_covs_present = [c for c in CAT_COVS if c in adata_hvg.obs.columns]
if len(cat_covs_present) < len(CAT_COVS):
    log.warning("  Missing covariate columns (skipped): %s",
                set(CAT_COVS) - set(cat_covs_present))

scvi.model.SCVI.setup_anndata(
    adata_hvg,
    layer                      = "counts",
    batch_key                  = BATCH_KEY,
    categorical_covariate_keys = cat_covs_present or None,
)


###############################################################################
# STEP 5: Train scVI (restart-safe via model checkpoint)
###############################################################################

_model_pt = os.path.join(SCVI_MODEL, "model.pt")
if os.path.exists(_model_pt):
    log.info("scVI model checkpoint found — loading from %s", SCVI_MODEL)
    model = scvi.model.SCVI.load(SCVI_MODEL, adata=adata_hvg)
else:
    model = scvi.model.SCVI(
        adata_hvg,
        n_latent        = N_LATENT,
        n_layers        = N_LAYERS,
        n_hidden        = N_HIDDEN,
        gene_likelihood = "nb",
        dispersion      = "gene-batch",
    )
    log.info("Model:\n%s", model)

    model.train(
        max_epochs              = N_EPOCHS,
        batch_size              = BATCH_SIZE,
        early_stopping          = True,
        early_stopping_patience = 20,
        plan_kwargs             = {"lr": LR},
        accelerator             = "gpu" if device == "cuda" else "cpu",
    )

    hist       = model.history
    elbo_train = hist["elbo_train"]
    elbo_val   = hist.get("elbo_validation", None)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(elbo_train.index, elbo_train["elbo_train"], label="train ELBO")
    if elbo_val is not None:
        ax.plot(elbo_val.index, elbo_val["elbo_validation"], label="val ELBO")
    ax.set_xlabel("Epoch"); ax.set_ylabel("ELBO"); ax.legend()
    ax.set_title("scVI training ELBO — Bladder Urothelium")
    fig.savefig(os.path.join(PLOT_DIR, "scvi_elbo.pdf"), bbox_inches="tight")
    plt.close(fig)

    model.save(SCVI_MODEL, overwrite=True)
    log.info("scVI model saved → %s", SCVI_MODEL)


###############################################################################
# STEP 6: Latent representation + UMAP + Leiden
###############################################################################

log.info("Extracting scVI latent representation ...")
adata_hvg.obsm["X_scVI"] = model.get_latent_representation()

log.info("Building neighbourhood graph ...")
sc.pp.neighbors(adata_hvg, use_rep="X_scVI", n_neighbors=N_NEIGHBOURS,
                key_added="scVI")

log.info("Leiden clustering (resolution=%.2f) ...", LEIDEN_RES)
sc.tl.leiden(adata_hvg, neighbors_key="scVI", resolution=LEIDEN_RES,
             key_added="leiden_scVI")

log.info("UMAP ...")
sc.tl.umap(adata_hvg, neighbors_key="scVI")
adata_hvg.obsm["X_umap_scVI"] = adata_hvg.obsm["X_umap"].copy()


###############################################################################
# STEP 7: scANVI (semi-supervised, optional)
###############################################################################

run_scanvi = False

if CELL_TYPE_KEY in adata_hvg.obs.columns:
    col = adata_hvg.obs[CELL_TYPE_KEY]
    if hasattr(col, "cat") and UNLABELED_CAT not in col.cat.categories:
        col = col.cat.add_categories(UNLABELED_CAT)
    adata_hvg.obs[CELL_TYPE_KEY] = col.fillna(UNLABELED_CAT).astype(str)
    n_labelled    = (adata_hvg.obs[CELL_TYPE_KEY] != UNLABELED_CAT).sum()
    frac_labelled = n_labelled / adata_hvg.n_obs
    log.info("Labelled cells for scANVI: %d / %d (%.1f %%)",
             n_labelled, adata_hvg.n_obs, 100 * frac_labelled)
    run_scanvi = frac_labelled >= SCANVI_LABEL_MIN_FRAC
    if not run_scanvi:
        log.info("  Too few labelled cells — skipping scANVI")
else:
    log.info("No '%s' column — skipping scANVI", CELL_TYPE_KEY)

if run_scanvi:
    log.info("Fine-tuning scANVI ...")
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        model,
        labels_key         = CELL_TYPE_KEY,
        unlabeled_category = UNLABELED_CAT,
    )
    scanvi_model.train(
        max_epochs          = 20,
        n_samples_per_label = 100,
        accelerator         = "gpu" if device == "cuda" else "cpu",
    )
    scanvi_model.save(SCANVI_MODEL, overwrite=True)
    log.info("scANVI model saved → %s", SCANVI_MODEL)

    adata_hvg.obsm["X_scANVI"]        = scanvi_model.get_latent_representation()
    adata_hvg.obs["cell_type_scanvi"] = scanvi_model.predict()

    sc.pp.neighbors(adata_hvg, use_rep="X_scANVI", n_neighbors=N_NEIGHBOURS,
                    key_added="scANVI")
    sc.tl.leiden(adata_hvg, neighbors_key="scANVI", resolution=LEIDEN_RES,
                 key_added="leiden_scANVI")
    sc.tl.umap(adata_hvg, neighbors_key="scANVI")
    adata_hvg.obsm["X_umap_scANVI"] = adata_hvg.obsm["X_umap"].copy()


###############################################################################
# STEP 8: UMAP plots
###############################################################################

def save_umap(adata, basis, color_keys, filename, title=""):
    present = [c for c in color_keys if c in adata.obs.columns]
    if not present:
        return
    sc.pl.embedding(adata, basis=basis, color=present,
                    ncols=2, show=False, return_fig=True)
    fig = plt.gcf()
    if title:
        fig.suptitle(title, y=1.01)
    fig.savefig(os.path.join(PLOT_DIR, filename), bbox_inches="tight", dpi=150)
    plt.close(fig)

umap_colors = [BATCH_KEY, "sample_id", "technology", "condition", "leiden_scVI"]
log.info("Saving UMAP plots ...")
save_umap(adata_hvg, "X_umap_scVI", umap_colors,
          "umap_scVI.pdf", "scVI — Bladder Urothelium")

if run_scanvi:
    save_umap(adata_hvg, "X_umap_scANVI",
              umap_colors + ["cell_type_scanvi"],
              "umap_scANVI.pdf", "scANVI — Bladder Urothelium")


###############################################################################
# STEP 9: Save
###############################################################################

adata_hvg.X = adata_hvg.layers["log_norm"]

log.info("Saving → %s", H5AD_OUT)
adata_hvg.write_h5ad(H5AD_OUT)

log.info("\n=== scVI integration complete — Bladder Urothelium ===")
log.info("  Cells          : %d", adata_hvg.n_obs)
log.info("  Genes (HVG)    : %d", adata_hvg.n_vars)
log.info("  scVI clusters  : %d", adata_hvg.obs["leiden_scVI"].nunique())
if run_scanvi:
    log.info("  scANVI clusters: %d",
             adata_hvg.obs["leiden_scANVI"].nunique())
log.info("  Output         : %s", H5AD_OUT)
