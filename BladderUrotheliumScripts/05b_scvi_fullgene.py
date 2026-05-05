#!/usr/bin/env python3
"""
05b_scvi_fullgene.py — Full-gene scVI integration for BladderUrothelium.

Strategy (mirrors RenalUrotheliumScripts steps 10+11):
  1. Load all per-sample qc_h5ad files; keep ALL genes (no HVG subset saved).
  2. Concatenate into one full-gene AnnData (outer join, fill NaN → 0).
  3. Normalize + log1p → layers["lognorm"].
  4. Select 3 000 batch-aware HVGs; train scVI on HVG subset.
  5. Get latent representation for ALL cells from HVG subset.
  6. Attach X_scVI embeddings to the FULL-gene object.
  7. Build neighbourhood graph, compute UMAP + Leiden.
  8. Save full-gene AnnData (all genes + scVI embeddings).

Input  : ../qc_h5ad/*.h5ad  (counts in X, per-sample QC'd)
Output : output/BladderUrothelium_allcells_scvi.h5ad  (~25 k genes)
         output/bladder_fullgene_scvi_model/
         output/bladder_fullgene_scvi_plots/
"""

import gc
import glob
import logging
import os

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
DATA_DIR   = "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/BladderUrothelium"
SCRIPT_DIR = os.path.join(DATA_DIR, "BladderUrotheliumScripts")
QC_DIR     = os.path.join(DATA_DIR, "qc_h5ad")
OUT_DIR    = os.path.join(SCRIPT_DIR, "output")
PLOT_DIR   = os.path.join(OUT_DIR, "bladder_fullgene_scvi_plots")
MODEL_DIR  = os.path.join(OUT_DIR, "bladder_fullgene_scvi_model")
H5AD_OUT   = os.path.join(OUT_DIR, "BladderUrothelium_allcells_scvi.h5ad")

for d in (OUT_DIR, PLOT_DIR):
    os.makedirs(d, exist_ok=True)

os.environ.setdefault("MPLCONFIGDIR", f"/tmp/matplotlib-{os.getuid()}")
os.environ.setdefault("NUMBA_CACHE_DIR", f"/tmp/numba-{os.getuid()}")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)
os.makedirs(os.environ["NUMBA_CACHE_DIR"], exist_ok=True)

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

BATCH_KEY = "gsm_id"
CAT_COVS  = ["technology", "sample_id"]

# ── Logging ────────────────────────────────────────────────────────────────────
logging.basicConfig(
    level   = logging.INFO,
    format  = "%(asctime)s  %(levelname)s  %(message)s",
    datefmt = "%H:%M:%S",
)
log = logging.getLogger(__name__)
scvi.settings.verbosity = 2
scvi.settings.dl_num_workers = int(os.environ.get("SLURM_CPUS_PER_TASK", 4))

device = "cuda" if torch.cuda.is_available() else "cpu"
log.info("PyTorch device: %s", device)
if device == "cuda":
    log.info("  GPU: %s", torch.cuda.get_device_name(0))


###############################################################################
# STEP 1: Load per-sample h5ad files (full gene set)
###############################################################################

h5ad_files = sorted(glob.glob(os.path.join(QC_DIR, "*.h5ad")))
log.info("Found %d h5ad files in %s", len(h5ad_files), QC_DIR)

adatas  = []
skipped = []

for path in h5ad_files:
    sid = os.path.splitext(os.path.basename(path))[0]
    a   = sc.read_h5ad(path)

    # Promote integer X → counts layer
    if "counts" not in a.layers:
        _x    = a.X
        _data = _x.data if sp.issparse(_x) else _x.ravel()
        if len(_data) > 0 and np.all(_data == np.floor(_data)):
            a.layers["counts"] = a.X.copy()
        else:
            log.info("  [skip – no counts] %s", sid)
            skipped.append((sid, "no counts layer"))
            continue

    if BATCH_KEY not in a.obs.columns or a.obs[BATCH_KEY].isna().all():
        log.warning("  '%s' missing for %s — using sample name", BATCH_KEY, sid)
        a.obs[BATCH_KEY] = sid

    log.info("  [loaded] %-28s  %6d cells × %6d genes", sid, a.n_obs, a.n_vars)
    adatas.append(a)

log.info("Samples loaded: %d  skipped: %d", len(adatas), len(skipped))
if not adatas:
    raise RuntimeError("No samples with raw counts found in %s" % QC_DIR)


###############################################################################
# STEP 2: Concatenate (outer join; fill NaN → 0)
###############################################################################

log.info("Concatenating (outer join, all genes) ...")
adata = ad.concat(adatas, join="outer", merge="same", index_unique="-")
del adatas
gc.collect()

for layer in list(adata.layers.keys()):
    if sp.issparse(adata.layers[layer]):
        adata.layers[layer].data = np.nan_to_num(adata.layers[layer].data)
    else:
        adata.layers[layer] = np.nan_to_num(adata.layers[layer])

log.info("Full-gene object: %d cells × %d genes", adata.n_obs, adata.n_vars)


###############################################################################
# STEP 3: Normalize + log1p (full-gene object)
###############################################################################

log.info("Normalizing for HVG selection + downstream use ...")
adata.X = adata.layers["counts"].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers["lognorm"] = adata.X.copy()


###############################################################################
# STEP 4: Select HVGs + build scVI subset
###############################################################################

log.info("Selecting %d batch-aware HVGs (batch_key=%s) ...", N_HVG, BATCH_KEY)
sc.pp.highly_variable_genes(
    adata,
    n_top_genes = N_HVG,
    batch_key   = BATCH_KEY,
    flavor      = "seurat_v3",
    layer       = "counts",
    subset      = False,
)
log.info("  HVGs selected: %d", adata.var["highly_variable"].sum())

# HVG subset for scVI training only
adata_hvg     = adata[:, adata.var["highly_variable"]].copy()
adata_hvg.X   = adata_hvg.layers["counts"]
log.info("HVG subset: %d cells × %d genes", adata_hvg.n_obs, adata_hvg.n_vars)


###############################################################################
# STEP 5: Train scVI (or load from checkpoint)
###############################################################################

cat_covs_present = [c for c in CAT_COVS
                    if c in adata_hvg.obs.columns and adata_hvg.obs[c].nunique() > 1]

scvi.model.SCVI.setup_anndata(
    adata_hvg,
    layer                      = "counts",
    batch_key                  = BATCH_KEY,
    categorical_covariate_keys = cat_covs_present or None,
)

_model_pt = os.path.join(MODEL_DIR, "model.pt")
if os.path.exists(_model_pt):
    log.info("Checkpoint found — loading scVI model from %s", MODEL_DIR)
    model = scvi.model.SCVI.load(MODEL_DIR, adata=adata_hvg)
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

    # ELBO curve
    hist       = model.history
    elbo_train = hist["elbo_train"]
    elbo_val   = hist.get("elbo_validation", None)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(elbo_train.index, elbo_train["elbo_train"], label="train ELBO")
    if elbo_val is not None:
        ax.plot(elbo_val.index, elbo_val["elbo_validation"], label="val ELBO")
    ax.set_xlabel("Epoch"); ax.set_ylabel("ELBO"); ax.legend()
    ax.set_title("scVI training ELBO — BladderUrothelium full-gene")
    fig.savefig(os.path.join(PLOT_DIR, "scvi_elbo.pdf"), bbox_inches="tight")
    plt.close(fig)

    model.save(MODEL_DIR, overwrite=True)
    log.info("scVI model saved → %s", MODEL_DIR)


###############################################################################
# STEP 6: Attach scVI embeddings to full-gene object
###############################################################################

log.info("Computing scVI latent representation (%d cells) ...", adata.n_obs)
adata.obsm["X_scVI"] = model.get_latent_representation()
log.info("  X_scVI shape: %s", adata.obsm["X_scVI"].shape)

del model, adata_hvg
gc.collect()


###############################################################################
# STEP 7: Neighbourhood graph → UMAP + Leiden (from scVI latent)
###############################################################################

log.info("Building neighbourhood graph (n_neighbors=%d) ...", N_NEIGHBOURS)
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=N_NEIGHBOURS, key_added="scVI")

log.info("Leiden clustering (resolution=%.2f) ...", LEIDEN_RES)
sc.tl.leiden(adata, neighbors_key="scVI", resolution=LEIDEN_RES, key_added="leiden_scVI")
log.info("  Clusters: %d", adata.obs["leiden_scVI"].nunique())

log.info("Computing UMAP ...")
sc.tl.umap(adata, neighbors_key="scVI")
adata.obsm["X_umap_scVI"] = adata.obsm["X_umap"].copy()


###############################################################################
# STEP 8: UMAP plots
###############################################################################

color_cols = [c for c in [BATCH_KEY, "sample_id", "technology", "condition", "leiden_scVI"]
              if c in adata.obs.columns]
if color_cols:
    log.info("Saving UMAP plots ...")
    fig = sc.pl.embedding(adata, basis="X_umap_scVI", color=color_cols,
                          ncols=2, show=False, return_fig=True)
    fig.savefig(os.path.join(PLOT_DIR, "umap_scVI_fullgene.pdf"),
                bbox_inches="tight", dpi=150)
    plt.close(fig)


###############################################################################
# STEP 9: Save full-gene object
###############################################################################

# Set X to lognorm for storage
adata.X = adata.layers["lognorm"]

log.info("Saving → %s", H5AD_OUT)
adata.write_h5ad(H5AD_OUT)

log.info("\n===== BladderUrothelium full-gene scVI complete =====")
log.info("  Cells          : %d", adata.n_obs)
log.info("  Genes (total)  : %d", adata.n_vars)
log.info("  scVI clusters  : %d", adata.obs["leiden_scVI"].nunique())
log.info("  obsm keys      : %s", list(adata.obsm.keys()))
log.info("  Output         : %s", H5AD_OUT)
