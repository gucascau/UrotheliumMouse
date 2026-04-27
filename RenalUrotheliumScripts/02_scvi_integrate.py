#!/usr/bin/env python3
"""
02_scvi_integrate.py — Train scVI on the preprocessed concatenated h5ad,
optionally fine-tune scANVI if labelled cells are present, compute UMAP
and Leiden clusters, and save the integrated object.

Input  : output/RenalUrothelium_preprocessed.h5ad
Outputs: output/RenalUrothelium_integrated.h5ad
         output/scvi_model/
         output/scanvi_model/   (if scANVI runs)
         output/plots/          (ELBO + UMAP PDFs)
"""

import os
import logging
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR     = "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
SCRIPT_DIR   = os.path.join(BASE_DIR, "RenalUrotheliumScripts")
OUT_DIR      = os.path.join(SCRIPT_DIR, "output")
PLOT_DIR     = os.path.join(OUT_DIR, "plots")
SCVI_MODEL   = os.path.join(OUT_DIR, "scvi_model")
SCANVI_MODEL = os.path.join(OUT_DIR, "scanvi_model")

IN_PATH  = os.path.join(OUT_DIR, "RenalUrothelium_preprocessed.h5ad")
OUT_PATH = os.path.join(OUT_DIR, "RenalUrothelium_integrated.h5ad")

os.makedirs(PLOT_DIR, exist_ok=True)

# ── Parameters ────────────────────────────────────────────────────────────────
BATCH_KEY    = "sample_id"
CAT_COVS     = ["technology"]

N_LATENT     = 20
N_LAYERS     = 2
N_HIDDEN     = 256
N_EPOCHS     = 400
BATCH_SIZE   = 256
LR           = 1e-3
LEIDEN_RES   = 0.5
N_NEIGHBOURS = 15

CELL_TYPE_KEY         = "cell_type_original"
UNLABELED_CAT         = "Unknown"
SCANVI_LABEL_MIN_FRAC = 0.01

# ── Logging ───────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)
scvi.settings.verbosity = 2

# ── Device ────────────────────────────────────────────────────────────────────
device = "cuda" if torch.cuda.is_available() else "cpu"
log.info("PyTorch device: %s", device)
if device == "cuda":
    log.info("  GPU: %s", torch.cuda.get_device_name(0))
scvi.settings.dl_num_workers = int(os.environ.get("SLURM_CPUS_PER_TASK", 4))


def _coerce_is_primary_data(series):
    """Convert mixed CellxGene values into a nullable boolean column."""
    mapping = {
        True: True,
        False: False,
        1: True,
        0: False,
        "1": True,
        "0": False,
        "true": True,
        "false": False,
        "True": True,
        "False": False,
        "": pd.NA,
        "nan": pd.NA,
        "None": pd.NA,
        "none": pd.NA,
    }
    coerced = series.map(lambda x: mapping.get(x, pd.NA))
    return coerced.astype("boolean")


def sanitise_dataframe_for_h5ad(df, axis_name):
    """Rename reserved columns and coerce object metadata into HDF5-safe dtypes."""
    df = df.copy()

    if "_index" in df.columns:
        replacement = f"{axis_name}_index_original"
        suffix = 1
        while replacement in df.columns:
            replacement = f"{axis_name}_index_original_{suffix}"
            suffix += 1
        df = df.rename(columns={"_index": replacement})
        log.info("Renamed reserved column '_index' → '%s' in %s", replacement, axis_name)

    for col in df.columns:
        series = df[col]
        if col == "is_primary_data":
            df[col] = _coerce_is_primary_data(series)
            continue
        if series.dtype == "object":
            non_na = series.dropna()
            if non_na.empty:
                df[col] = series.fillna("").astype(str)
                continue

            if non_na.map(lambda x: isinstance(x, (bool, np.bool_))).all():
                df[col] = series.astype("boolean")
                continue

            if non_na.map(lambda x: isinstance(x, (int, float, np.integer, np.floating, bool, np.bool_))).all():
                df[col] = pd.to_numeric(series, errors="coerce")
                continue

            df[col] = series.fillna("").map(str)

    return df


def sanitise_anndata_for_h5ad(adata):
    """Normalise obs/var metadata so AnnData can be written reliably."""
    adata.obs = sanitise_dataframe_for_h5ad(adata.obs, "obs")
    adata.var = sanitise_dataframe_for_h5ad(adata.var, "var")
    return adata


###############################################################################
# STEP 1: Load preprocessed data
###############################################################################

log.info("Loading %s ...", IN_PATH)
adata = sc.read_h5ad(IN_PATH)
log.info("  %d cells × %d genes", adata.n_obs, adata.n_vars)
log.info("  Layers: %s", list(adata.layers.keys()))

# Preserve log-norm X for post-integration visualisation
adata.layers["log_norm"] = adata.X.copy()

# Set X to raw counts for scVI
if "counts" not in adata.layers:
    log.warning("No 'counts' layer — using X as counts (suboptimal).")
    adata.layers["counts"] = adata.X.copy()
adata.X = adata.layers["counts"]

log.info("  Samples: %s", sorted(adata.obs[BATCH_KEY].unique()))


###############################################################################
# STEP 2: scVI setup + training
###############################################################################

cat_covs = [c for c in CAT_COVS if c in adata.obs.columns
            and adata.obs[c].nunique() > 1]
if len(cat_covs) < len(CAT_COVS):
    log.warning("Skipped categorical covariates (missing or single-level): %s",
                set(CAT_COVS) - set(cat_covs))

log.info("Setting up scVI (batch=%s, covariates=%s) ...", BATCH_KEY, cat_covs)
scvi.model.SCVI.setup_anndata(
    adata,
    layer                      = "counts",
    batch_key                  = BATCH_KEY,
    categorical_covariate_keys = cat_covs or None,
)

model = scvi.model.SCVI(
    adata,
    n_latent        = N_LATENT,
    n_layers        = N_LAYERS,
    n_hidden        = N_HIDDEN,
    gene_likelihood = "nb",
    dispersion      = "gene-batch",
)
log.info("  Model summary:\n%s", model)

log.info("Training scVI (max_epochs=%d) ...", N_EPOCHS)
model.train(
    max_epochs              = N_EPOCHS,
    batch_size              = BATCH_SIZE,
    early_stopping          = True,
    early_stopping_patience = 20,
    plan_kwargs             = {"lr": LR},
    accelerator             = "gpu" if device == "cuda" else "cpu",
)

# Save ELBO plot
hist = model.history
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(hist["elbo_train"].index,
        hist["elbo_train"]["elbo_train"], label="train")
if "elbo_validation" in hist:
    ax.plot(hist["elbo_validation"].index,
            hist["elbo_validation"]["elbo_validation"], label="val")
ax.set_xlabel("Epoch"); ax.set_ylabel("ELBO"); ax.legend()
ax.set_title("scVI training ELBO — RenalUrothelium")
fig.savefig(os.path.join(PLOT_DIR, "scvi_elbo.pdf"), bbox_inches="tight")
plt.close(fig)

model.save(SCVI_MODEL, overwrite=True)
log.info("scVI model saved → %s", SCVI_MODEL)


###############################################################################
# STEP 3: Latent representation, neighbours, UMAP, Leiden
###############################################################################

log.info("Extracting scVI latent embedding ...")
adata.obsm["X_scVI"] = model.get_latent_representation()

log.info("Building neighbourhood graph ...")
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=N_NEIGHBOURS,
                key_added="scVI")

log.info("Leiden clustering (resolution=%.2f) ...", LEIDEN_RES)
sc.tl.leiden(adata, neighbors_key="scVI", resolution=LEIDEN_RES,
             key_added="leiden_scVI")

log.info("UMAP on scVI latent ...")
sc.tl.umap(adata, neighbors_key="scVI")
adata.obsm["X_umap_scVI"] = adata.obsm["X_umap"].copy()


###############################################################################
# STEP 4: scANVI (semi-supervised, optional)
###############################################################################

run_scanvi = False

if CELL_TYPE_KEY in adata.obs.columns:
    adata.obs[CELL_TYPE_KEY] = (
        adata.obs[CELL_TYPE_KEY].fillna(UNLABELED_CAT).astype(str)
    )
    n_lab  = (adata.obs[CELL_TYPE_KEY] != UNLABELED_CAT).sum()
    frac   = n_lab / adata.n_obs
    log.info("Labelled cells: %d / %d (%.1f%%)", n_lab, adata.n_obs, 100*frac)
    run_scanvi = frac >= SCANVI_LABEL_MIN_FRAC
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

    adata.obsm["X_scANVI"]        = scanvi_model.get_latent_representation()
    adata.obs["cell_type_scanvi"] = scanvi_model.predict()

    sc.pp.neighbors(adata, use_rep="X_scANVI", n_neighbors=N_NEIGHBOURS,
                    key_added="scANVI")
    sc.tl.leiden(adata, neighbors_key="scANVI", resolution=LEIDEN_RES,
                 key_added="leiden_scANVI")
    sc.tl.umap(adata, neighbors_key="scANVI")
    adata.obsm["X_umap_scANVI"] = adata.obsm["X_umap"].copy()
    log.info("scANVI done → %s", SCANVI_MODEL)


###############################################################################
# STEP 5: UMAP plots
###############################################################################

def plot_umap(embedding, colors, fname, title=""):
    cols = [c for c in colors if c in adata.obs.columns]
    if not cols:
        return
    fig = sc.pl.embedding(adata, basis=embedding, color=cols,
                          ncols=2, show=False, return_fig=True)
    if title:
        fig.suptitle(title, y=1.01)
    fig.savefig(os.path.join(PLOT_DIR, fname), bbox_inches="tight", dpi=150)
    plt.close(fig)

base_colors = [BATCH_KEY, "technology", "condition", "leiden_scVI"]
plot_umap("X_umap_scVI", base_colors, "umap_scVI.pdf", "scVI — RenalUrothelium")

if run_scanvi:
    plot_umap("X_umap_scANVI",
              base_colors + ["cell_type_scanvi", "leiden_scANVI"],
              "umap_scANVI.pdf", "scANVI — RenalUrothelium")


###############################################################################
# STEP 6: Save
###############################################################################

# Restore log-norm X for downstream scanpy / R use
adata.X = adata.layers["log_norm"]
adata = sanitise_anndata_for_h5ad(adata)

log.info("Saving integrated AnnData → %s", OUT_PATH)
adata.write_h5ad(OUT_PATH)

log.info("\n===== scVI integration complete =====")
log.info("  Cells         : %d", adata.n_obs)
log.info("  scVI clusters : %d", adata.obs["leiden_scVI"].nunique())
if run_scanvi:
    log.info("  scANVI clusters: %d", adata.obs["leiden_scANVI"].nunique())
log.info("  Output        : %s", OUT_PATH)
