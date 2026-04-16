#!/usr/bin/env python3
"""
Script 06: scVI / scANVI integration from per-sample QC h5ad files

Workflow
--------
1. Load all *.h5ad from qc_h5ad/ and retain only samples that carry a
   raw 'counts' layer (scVI requires integer counts; pre-normalised-only
   samples such as ChenSpatial, LakesnRNA, MKA, KudoUUOUrothelium,
   KidneyUUO7/8 are skipped automatically).
2. Concatenate into a single AnnData (outer join on genes; missing values
   filled with 0).
3. Exclude non-kidney samples (BladderHomogenate1/2).
4. Normalize + log1p on X for HVG selection; keep raw counts in layers.
5. Select HVGs with scanpy (batch-aware, flavour="seurat_v3").
6. Subset to HVGs and train scVI (batch_key = sample_id,
   categorical covariate = technology).
7. Extract latent Z, build neighbours, run UMAP + Leiden clustering.
8. Optionally fine-tune scANVI if cell_type_original labels are available.
9. Save integrated AnnData + model directories.

Outputs (under integration_output/)
-------------------------------------
  scvi_qc_integrated.h5ad  — integrated AnnData (X = log-norm, obsm = embeddings)
  scvi_qc_model/           — trained scVI model
  scanvi_qc_model/         — trained scANVI model (if run)
  scvi_qc_plots/           — UMAP PDFs + ELBO curve

Requirements:
  scvi-tools>=1.1  scanpy>=1.9  anndata  torch  leidenalg  igraph  matplotlib
"""

import os
import glob
import logging
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scvi
import torch
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR     = "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
QC_H5AD_DIR  = os.path.join(DATA_DIR, "qc_h5ad")
OUT_DIR      = os.path.join(DATA_DIR, "integration_output")
PLOT_DIR     = os.path.join(OUT_DIR, "scvi_qc_plots")
SCVI_MODEL   = os.path.join(OUT_DIR, "scvi_qc_model")
SCANVI_MODEL = os.path.join(OUT_DIR, "scanvi_qc_model")
H5AD_OUT     = os.path.join(OUT_DIR, "scvi_qc_integrated.h5ad")

for d in (OUT_DIR, PLOT_DIR):
    os.makedirs(d, exist_ok=True)

# ── Parameters ─────────────────────────────────────────────────────────────────
N_HVG        = 3000    # number of highly variable genes
N_LATENT     = 20      # scVI latent dimensions
N_LAYERS     = 2       # encoder/decoder depth
N_HIDDEN     = 256     # neurons per hidden layer
N_EPOCHS     = 400     # max training epochs (early stopping active)
BATCH_SIZE   = 256
LR           = 1e-3
LEIDEN_RES   = 0.5
N_NEIGHBOURS = 15

BATCH_KEY    = "gsm_id"       # correct for per-GSM/GSE technical batch effects
CAT_COVS     = ["technology", "sample_id"]  # additional categorical covariates

CELL_TYPE_KEY         = "cell_type_original"
UNLABELED_CAT         = "Unknown"
SCANVI_LABEL_MIN_FRAC = 0.01

# Samples to exclude from integration (non-kidney tissue)
EXCLUDE_SAMPLES = {"BladderHomogenate1", "BladderHomogenate2"}

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

adatas      = []
skipped     = []

for path in h5ad_files:
    sid = os.path.splitext(os.path.basename(path))[0]

    # Skip non-kidney samples
    if sid in EXCLUDE_SAMPLES:
        log.info("  [skip – excluded] %s", sid)
        skipped.append((sid, "excluded"))
        continue

    a = sc.read_h5ad(path)

    # scVI requires raw integer counts — skip samples that only have
    # a normalised data layer (ChenSpatial, LakesnRNA, MKA, KudoUUOUrothelium,
    # KidneyUUO7, KidneyUUO8).
    if "counts" not in a.layers:
        log.info("  [skip – no counts] %s  (%d cells)", sid, a.n_obs)
        skipped.append((sid, "no counts layer"))
        continue

    # Ensure gsm_id is present; fall back to sample_id for datasets that
    # were not deposited under a GSM accession.
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


###############################################################################
# STEP 2: Concatenate
###############################################################################

log.info("\nConcatenating …")
# join="outer" fills missing genes with 0 so all samples share a common
# gene universe after concatenation.
adata = ad.concat(
    adatas,
    join        = "outer",
    merge       = "same",   # keep obs/var columns that are identical across all
    label       = BATCH_KEY,
    keys        = [a.obs[BATCH_KEY].iloc[0] for a in adatas],
    index_unique= "-",
)

# Fill NaN counts (from outer join on missing genes) with 0
import scipy.sparse as sp
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

# Store raw counts before normalization
adata.layers["counts"] = adata.layers["counts"].copy()

# Normalize X for HVG selection (counts → log-norm)
log.info("Normalizing for HVG selection …")
adata.X = adata.layers["counts"].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers["log_norm"] = adata.X.copy()

# Batch-aware HVG selection using seurat_v3 (works on raw counts)
log.info("Selecting %d HVGs (batch_key=%s) …", N_HVG, BATCH_KEY)
sc.pp.highly_variable_genes(
    adata,
    n_top_genes = N_HVG,
    batch_key   = BATCH_KEY,
    flavor      = "seurat_v3",
    layer       = "counts",
    subset      = False,   # keep all genes; flag only
)
log.info("  HVGs selected: %d", adata.var["highly_variable"].sum())


###############################################################################
# STEP 4: Subset to HVGs + set up scVI
###############################################################################

adata_hvg = adata[:, adata.var["highly_variable"]].copy()
log.info("After HVG subset: %d cells × %d genes", adata_hvg.n_obs, adata_hvg.n_vars)

# X must be raw counts for scVI
adata_hvg.X = adata_hvg.layers["counts"]

log.info("Setting up scVI …")
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
# STEP 5: Train scVI
###############################################################################

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

# ELBO training curve
hist       = model.history
elbo_train = hist["elbo_train"]
elbo_val   = hist.get("elbo_validation", None)
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(elbo_train.index, elbo_train["elbo_train"], label="train ELBO")
if elbo_val is not None:
    ax.plot(elbo_val.index, elbo_val["elbo_validation"], label="val ELBO")
ax.set_xlabel("Epoch"); ax.set_ylabel("ELBO"); ax.legend()
ax.set_title("scVI training ELBO")
fig.savefig(os.path.join(PLOT_DIR, "scvi_elbo.pdf"), bbox_inches="tight")
plt.close(fig)

model.save(SCVI_MODEL, overwrite=True)
log.info("scVI model saved → %s", SCVI_MODEL)


###############################################################################
# STEP 6: Latent representation + UMAP + Leiden
###############################################################################

log.info("Extracting scVI latent representation …")
adata_hvg.obsm["X_scVI"] = model.get_latent_representation()

log.info("Building neighbourhood graph …")
sc.pp.neighbors(adata_hvg, use_rep="X_scVI", n_neighbors=N_NEIGHBOURS,
                key_added="scVI")

log.info("Leiden clustering (resolution=%.2f) …", LEIDEN_RES)
sc.tl.leiden(adata_hvg, neighbors_key="scVI", resolution=LEIDEN_RES,
             key_added="leiden_scVI")

log.info("UMAP …")
sc.tl.umap(adata_hvg, neighbors_key="scVI")
adata_hvg.obsm["X_umap_scVI"] = adata_hvg.obsm["X_umap"].copy()


###############################################################################
# STEP 7: scANVI (semi-supervised, optional)
###############################################################################

run_scanvi = False

if CELL_TYPE_KEY is not None and CELL_TYPE_KEY in adata_hvg.obs.columns:
    adata_hvg.obs[CELL_TYPE_KEY] = (
        adata_hvg.obs[CELL_TYPE_KEY].fillna(UNLABELED_CAT).astype(str)
    )
    n_labelled    = (adata_hvg.obs[CELL_TYPE_KEY] != UNLABELED_CAT).sum()
    frac_labelled = n_labelled / adata_hvg.n_obs
    log.info("Labelled cells for scANVI: %d / %d (%.1f %%)",
             n_labelled, adata_hvg.n_obs, 100 * frac_labelled)

    if frac_labelled >= SCANVI_LABEL_MIN_FRAC:
        run_scanvi = True
    else:
        log.info("  Too few labelled cells — skipping scANVI")
else:
    log.info("No '%s' column — skipping scANVI", CELL_TYPE_KEY)

if run_scanvi:
    log.info("Fine-tuning scANVI …")
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
log.info("Saving UMAP plots …")
save_umap(adata_hvg, "X_umap_scVI", umap_colors,
          "umap_scVI.pdf", "scVI integration (QC samples)")

if run_scanvi:
    save_umap(adata_hvg, "X_umap_scANVI",
              umap_colors + ["cell_type_scanvi"],
              "umap_scANVI.pdf", "scANVI integration (QC samples)")


###############################################################################
# STEP 9: Save
###############################################################################

# Switch X back to log-normalised for downstream scanpy analysis
adata_hvg.X = adata_hvg.layers["log_norm"]

log.info("Saving → %s", H5AD_OUT)
adata_hvg.write_h5ad(H5AD_OUT)

log.info("\n=== scVI/scANVI integration complete ===")
log.info("  Cells          : %d", adata_hvg.n_obs)
log.info("  Genes (HVG)    : %d", adata_hvg.n_vars)
log.info("  scVI clusters  : %d", adata_hvg.obs["leiden_scVI"].nunique())
if run_scanvi:
    log.info("  scANVI clusters: %d",
             adata_hvg.obs["leiden_scANVI"].nunique())
log.info("  Output         : %s", H5AD_OUT)
