#!/usr/bin/env python3
"""
Script 04: scVI / scANVI integration of mouse kidney UUO single-cell data

Workflow
--------
1. Load h5ad exported by 03b_export_for_scvi.R
2. Subset to HVGs, log-normalise a copy for PCA/UMAP visualisation
3. Train scVI (VAE) with batch_key = "sample_id" and
   categorical covariate "technology" to correct cross-platform effects
4. Get latent representation Z (50-dim by default)
5. Run UMAP and Leiden clustering on Z
6. If cell_type_original labels exist for ≥1 % of cells, fine-tune
   scANVI (semi-supervised) and predict labels for unlabelled cells
7. Save:
     integration_output/scvi_integrated.h5ad  — full AnnData with all embeddings
     integration_output/scvi_model/            — trained scVI model directory
     integration_output/scanvi_model/          — trained scANVI model directory (if run)

Requirements (install via pip or conda):
  scvi-tools>=1.1  scanpy>=1.9  torch  leidenalg  igraph  matplotlib

GPU: strongly recommended for ≥200k cells.  The script auto-detects CUDA.
"""

import os
import sys
import logging
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch
import matplotlib
matplotlib.use("Agg")   # headless rendering on HPC
import matplotlib.pyplot as plt

# ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR = "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
OUT_DIR  = os.path.join(DATA_DIR, "integration_output")
os.makedirs(OUT_DIR, exist_ok=True)

H5AD_IN         = os.path.join(OUT_DIR, "scvi_input.h5ad")
H5AD_OUT        = os.path.join(OUT_DIR, "scvi_integrated.h5ad")
SCVI_MODEL_DIR  = os.path.join(OUT_DIR, "scvi_model")
SCANVI_MODEL_DIR= os.path.join(OUT_DIR, "scanvi_model")
PLOT_DIR        = os.path.join(OUT_DIR, "scvi_plots")
os.makedirs(PLOT_DIR, exist_ok=True)

# ── Parameters ─────────────────────────────────────────────────────────────────
N_LATENT     = 20      # scVI latent dimensions (10–30 typical for large datasets)
N_LAYERS     = 2       # encoder/decoder depth
N_HIDDEN     = 256     # neurons per hidden layer
N_EPOCHS     = 400     # max training epochs (early stopping active)
BATCH_SIZE   = 256     # increase to 512 on multi-GPU nodes
LR           = 1e-3    # initial learning rate
LEIDEN_RES   = 0.5     # Leiden clustering resolution on scVI neighbours
N_NEIGHBOURS = 15      # k for scVI neighbourhood graph

# Covariates
BATCH_KEY    = "sample_id"            # main integration batch
CAT_COVS     = ["technology"]         # additional categorical covariates
# CON_COVS   = ["pct_mt"]            # uncomment to regress continuous covariates

# scANVI: treat cells with this label string as unlabelled
UNLABELED_CAT = "Unknown"
CELL_TYPE_KEY = "cell_type_original"  # column in obs; set None to skip scANVI
# Minimum fraction of cells that must be labelled to trigger scANVI training
SCANVI_LABEL_MIN_FRAC = 0.01

# ── Logging ────────────────────────────────────────────────────────────────────
logging.basicConfig(
    level   = logging.INFO,
    format  = "%(asctime)s  %(levelname)s  %(message)s",
    datefmt = "%H:%M:%S",
)
log = logging.getLogger(__name__)
scvi.settings.verbosity = 2   # INFO level from scvi-tools

# ── Device ─────────────────────────────────────────────────────────────────────
device = "cuda" if torch.cuda.is_available() else "cpu"
log.info("PyTorch device: %s", device)
if device == "cuda":
    log.info("  GPU: %s", torch.cuda.get_device_name(0))
scvi.settings.dl_num_workers = int(os.environ.get("SLURM_CPUS_PER_TASK", 4))


###############################################################################
# STEP 1: Load data
###############################################################################

log.info("Loading %s …", H5AD_IN)
adata = sc.read_h5ad(H5AD_IN)
log.info("  Loaded: %d cells × %d genes", adata.n_obs, adata.n_vars)
log.info("  Samples: %s", sorted(adata.obs[BATCH_KEY].unique()))

# Store raw counts before any normalisation step
adata.layers["counts"] = adata.X.copy()


###############################################################################
# STEP 2: Subset to HVGs  +  log-normalise for PCA / UMAP visualisation
###############################################################################

log.info("Subsetting to %d HVGs …", adata.var["highly_variable"].sum())
adata = adata[:, adata.var["highly_variable"]].copy()
log.info("  After subsetting: %d cells × %d genes", adata.n_obs, adata.n_vars)

# Log-normalise a copy for visualisation (stored separately; scVI uses raw counts)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers["log_norm"] = adata.X.copy()

# Set X back to raw counts for scVI
adata.X = adata.layers["counts"]


###############################################################################
# STEP 3: scVI setup + training
###############################################################################

log.info("Setting up scVI …")

# Validate that CAT_COVS columns exist
cat_covs_present = [c for c in CAT_COVS if c in adata.obs.columns]
if len(cat_covs_present) < len(CAT_COVS):
    missing = set(CAT_COVS) - set(cat_covs_present)
    log.warning("  Missing covariate columns (skipped): %s", missing)

scvi.model.SCVI.setup_anndata(
    adata,
    layer                      = "counts",
    batch_key                  = BATCH_KEY,
    categorical_covariate_keys = cat_covs_present if cat_covs_present else None,
)

model = scvi.model.SCVI(
    adata,
    n_latent  = N_LATENT,
    n_layers  = N_LAYERS,
    n_hidden  = N_HIDDEN,
    gene_likelihood = "nb",          # negative binomial — appropriate for scRNA-seq counts
    dispersion      = "gene-batch",  # per-gene, per-batch overdispersion
)
log.info("  Model:\n%s", model)

model.train(
    max_epochs         = N_EPOCHS,
    batch_size         = BATCH_SIZE,
    early_stopping     = True,
    early_stopping_patience = 20,
    plan_kwargs        = {"lr": LR},
    accelerator        = "gpu" if device == "cuda" else "cpu",
)

# Save training ELBO history
train_hist = model.history
elbo_train = train_hist["elbo_train"]
elbo_val   = train_hist.get("elbo_validation", None)

fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(elbo_train.index, elbo_train["elbo_train"], label="train ELBO")
if elbo_val is not None:
    ax.plot(elbo_val.index, elbo_val["elbo_validation"], label="val ELBO")
ax.set_xlabel("Epoch"); ax.set_ylabel("ELBO"); ax.legend()
ax.set_title("scVI training ELBO")
fig.savefig(os.path.join(PLOT_DIR, "scvi_elbo.pdf"), bbox_inches="tight")
plt.close(fig)

model.save(SCVI_MODEL_DIR, overwrite=True)
log.info("scVI model saved → %s", SCVI_MODEL_DIR)


###############################################################################
# STEP 4: Latent representation + downstream analysis
###############################################################################

log.info("Extracting scVI latent representation …")
adata.obsm["X_scVI"] = model.get_latent_representation()

log.info("Building neighbourhood graph on X_scVI …")
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=N_NEIGHBOURS, key_added="scVI")

log.info("Leiden clustering (resolution=%.2f) …", LEIDEN_RES)
sc.tl.leiden(adata, neighbors_key="scVI", resolution=LEIDEN_RES,
             key_added="leiden_scVI")

log.info("UMAP on scVI latent …")
sc.tl.umap(adata, neighbors_key="scVI")
adata.obsm["X_umap_scVI"] = adata.obsm["X_umap"].copy()


###############################################################################
# STEP 5: scANVI (semi-supervised) — optional
###############################################################################

run_scanvi = False

if CELL_TYPE_KEY is not None and CELL_TYPE_KEY in adata.obs.columns:
    # Fill NaN labels with the unlabelled category string
    adata.obs[CELL_TYPE_KEY] = (
        adata.obs[CELL_TYPE_KEY]
        .fillna(UNLABELED_CAT)
        .astype(str)
    )
    n_labelled = (adata.obs[CELL_TYPE_KEY] != UNLABELED_CAT).sum()
    frac_labelled = n_labelled / adata.n_obs
    log.info("Labelled cells for scANVI: %d / %d (%.1f %%)",
             n_labelled, adata.n_obs, 100 * frac_labelled)

    if frac_labelled >= SCANVI_LABEL_MIN_FRAC:
        run_scanvi = True
    else:
        log.info("  Too few labelled cells — skipping scANVI")
else:
    log.info("No '%s' column found — skipping scANVI", CELL_TYPE_KEY)

if run_scanvi:
    log.info("Fine-tuning scANVI from scVI model …")
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        model,
        labels_key         = CELL_TYPE_KEY,
        unlabeled_category = UNLABELED_CAT,
    )
    scanvi_model.train(
        max_epochs        = 20,
        n_samples_per_label = 100,
        accelerator       = "gpu" if device == "cuda" else "cpu",
    )
    scanvi_model.save(SCANVI_MODEL_DIR, overwrite=True)
    log.info("scANVI model saved → %s", SCANVI_MODEL_DIR)

    adata.obsm["X_scANVI"]          = scanvi_model.get_latent_representation()
    adata.obs["cell_type_scanvi"]   = scanvi_model.predict()

    # UMAP on scANVI latent
    sc.pp.neighbors(adata, use_rep="X_scANVI", n_neighbors=N_NEIGHBOURS,
                    key_added="scANVI")
    sc.tl.leiden(adata, neighbors_key="scANVI", resolution=LEIDEN_RES,
                 key_added="leiden_scANVI")
    sc.tl.umap(adata, neighbors_key="scANVI")
    adata.obsm["X_umap_scANVI"] = adata.obsm["X_umap"].copy()


###############################################################################
# STEP 6: Plots
###############################################################################

def save_umap(adata, embedding_key, color_keys, filename, title=""):
    sc.pl.embedding(
        adata, basis=embedding_key, color=color_keys,
        ncols=2, show=False, return_fig=True
    )
    fig = plt.gcf()
    if title:
        fig.suptitle(title, y=1.01)
    fig.savefig(os.path.join(PLOT_DIR, filename), bbox_inches="tight", dpi=150)
    plt.close(fig)

umap_colors = [BATCH_KEY, "technology", "condition", "leiden_scVI"]
umap_colors = [c for c in umap_colors if c in adata.obs.columns]

log.info("Saving UMAP plots …")
save_umap(adata, "X_umap_scVI", umap_colors,
          "umap_scVI.pdf", "scVI integration")

if run_scanvi:
    scanvi_colors = umap_colors + ["cell_type_scanvi"]
    scanvi_colors = [c for c in scanvi_colors if c in adata.obs.columns]
    save_umap(adata, "X_umap_scANVI", scanvi_colors,
              "umap_scANVI.pdf", "scANVI integration")


###############################################################################
# STEP 7: Save AnnData
###############################################################################

# Switch X back to log-normalised values for downstream analysis in scanpy
adata.X = adata.layers["log_norm"]

log.info("Saving integrated AnnData → %s", H5AD_OUT)
adata.write_h5ad(H5AD_OUT)

log.info("\n=== scVI/scANVI integration complete ===")
log.info("  Cells          : %d", adata.n_obs)
log.info("  scVI clusters  : %d", adata.obs["leiden_scVI"].nunique())
if run_scanvi:
    log.info("  scANVI clusters: %d", adata.obs["leiden_scANVI"].nunique())
    log.info("  scANVI predicted types: %s",
             sorted(adata.obs["cell_type_scanvi"].unique()))
log.info("  Output         : %s", OUT_DIR)
