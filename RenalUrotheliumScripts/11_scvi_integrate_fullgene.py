#!/usr/bin/env python3
"""
11_scvi_integrate_fullgene.py — Attach scVI embeddings to the full-gene
object by reusing the already-trained scVI model (output/scvi_model/).

No re-training is performed.  The existing model was trained on the 3,000
batch-aware HVGs across all 1.1M cells.  This script:
  1. Loads RenalUrothelium_allcells_fullgene.h5ad  (~25k genes, 1.1M cells).
  2. Subsets to the same 3,000 HVGs the model knows about.
  3. Loads the existing scVI model and calls get_latent_representation().
  4. Computes UMAP and Leiden clusters from the latent space.
  5. Transfers scanvi_Lake/MKA annotations from the integrated annotated object.
  6. Saves RenalUrothelium_allcells_scvi.h5ad with all ~25k genes retained
     plus the scVI UMAP and cluster labels.

Usage
-----
  python 11_scvi_integrate_fullgene.py
"""

import logging
import os
import scanpy as sc
import scvi
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR   = "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
SCRIPT_DIR = os.path.join(BASE_DIR, "RenalUrotheliumScripts")
OUT_DIR    = os.path.join(SCRIPT_DIR, "output")
PLOT_DIR   = os.path.join(OUT_DIR, "plots")

IN_PATH    = os.path.join(OUT_DIR, "RenalUrothelium_allcells_fullgene.h5ad")
OUT_PATH   = os.path.join(OUT_DIR, "RenalUrothelium_allcells_scvi.h5ad")
MODEL_DIR   = os.path.join(OUT_DIR, "scvi_model")          # existing trained model
SCANVI_DIR  = os.path.join(OUT_DIR, "scanvi_model")        # existing scANVI model
ANNOT_H5AD = os.path.join(OUT_DIR, "RenalUrothelium_integrated_annotated.h5ad")

os.makedirs(PLOT_DIR, exist_ok=True)

# ── Integration parameters (must match original 02_scvi_integrate.py) ─────────
BATCH_KEY    = "sample_id"
CAT_COVS     = ["technology"]
N_NEIGHBOURS = 15
LEIDEN_RES   = 0.5

logging.basicConfig(
    level  = logging.INFO,
    format = "%(asctime)s  %(levelname)s  %(message)s",
    datefmt= "%H:%M:%S",
)
log = logging.getLogger(__name__)


def main():
    # ── Step 1: Load full-gene object ─────────────────────────────────────────
    log.info("Loading full-gene object: %s", IN_PATH)
    adata = sc.read_h5ad(IN_PATH)
    log.info("  %d cells × %d genes", adata.n_obs, adata.n_vars)

    # ── Step 2: Get the 3,000 HVG names the model was trained on ──────────────
    log.info("Reading HVG gene list from integrated annotated object ...")
    adata_int = sc.read_h5ad(ANNOT_H5AD, backed="r")
    hvg_genes = list(adata_int.var_names)          # the 3,000 HVGs
    log.info("  HVGs in existing model: %d", len(hvg_genes))
    adata_int.file.close()

    # Subset full-gene object to those HVGs for model loading
    hvg_present = [g for g in hvg_genes if g in adata.var_names]
    missing     = len(hvg_genes) - len(hvg_present)
    if missing > 0:
        log.warning("  %d HVGs not found in full-gene object (gene name mismatch?) — "
                    "using %d / %d", missing, len(hvg_present), len(hvg_genes))

    adata_hvg = adata[:, hvg_present].copy()
    log.info("  Subset for model: %d cells × %d genes", adata_hvg.n_obs, adata_hvg.n_vars)

    # ── Step 3: Load existing scVI model ──────────────────────────────────────
    log.info("Loading scVI model from %s ...", MODEL_DIR)

    cat_covs = [c for c in CAT_COVS if c in adata_hvg.obs.columns
                and adata_hvg.obs[c].nunique() > 1]

    # setup_anndata must match exactly what was used during training
    scvi.model.SCVI.setup_anndata(
        adata_hvg,
        layer                      = "counts",
        batch_key                  = BATCH_KEY,
        categorical_covariate_keys = cat_covs or None,
    )

    model = scvi.model.SCVI.load(MODEL_DIR, adata=adata_hvg)
    log.info("  Model loaded successfully")

    # ── Step 4a: scVI latent representation ───────────────────────────────────
    log.info("Computing scVI latent representation (%d cells) ...", adata.n_obs)
    adata.obsm["X_scVI"] = model.get_latent_representation()
    log.info("  X_scVI shape: %s", adata.obsm["X_scVI"].shape)

    # ── Step 4b: scANVI latent representation (if model exists) ───────────────
    if os.path.isdir(SCANVI_DIR):
        log.info("Loading scANVI model from %s ...", SCANVI_DIR)
        try:
            scanvi_model = scvi.model.SCANVI.load(SCANVI_DIR, adata=adata_hvg)
            log.info("  scANVI model loaded")
            adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation()
            log.info("  X_scANVI shape: %s", adata.obsm["X_scANVI"].shape)

            sc.pp.neighbors(adata, use_rep="X_scANVI", n_neighbors=N_NEIGHBOURS,
                            key_added="scANVI")
            sc.tl.leiden(adata, neighbors_key="scANVI", resolution=LEIDEN_RES,
                         key_added="leiden_scANVI")
            sc.tl.umap(adata, neighbors_key="scANVI")
            adata.obsm["X_umap_scANVI"] = adata.obsm["X_umap"].copy()
            log.info("  scANVI clusters: %d", adata.obs["leiden_scANVI"].nunique())
        except Exception as e:
            log.warning("  scANVI loading failed (%s) — skipping", e)
    else:
        log.info("No scANVI model found at %s — skipping", SCANVI_DIR)

    # ── Step 5: Neighbourhood graph, UMAP, Leiden (from scVI) ─────────────────
    log.info("Building neighbourhood graph from scVI latent ...")
    sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=N_NEIGHBOURS,
                    key_added="scVI")

    log.info("Leiden clustering (resolution=%.2f) ...", LEIDEN_RES)
    sc.tl.leiden(adata, neighbors_key="scVI", resolution=LEIDEN_RES,
                 key_added="leiden_scVI")

    log.info("Computing UMAP from scVI latent ...")
    sc.tl.umap(adata, neighbors_key="scVI")
    adata.obsm["X_umap_scVI"] = adata.obsm["X_umap"].copy()

    log.info("  scVI clusters: %d", adata.obs["leiden_scVI"].nunique())

    # ── Step 6: Transfer annotations from integrated annotated object ──────────
    annot_cols = ["scanvi_Lake_label", "scanvi_Lake_confidence",
                  "scanvi_MKA_label",  "scanvi_MKA_confidence",
                  "leiden_scVI"]
    log.info("Transferring annotations from %s ...", ANNOT_H5AD)
    adata_annot = sc.read_h5ad(ANNOT_H5AD, backed="r")
    available   = [c for c in annot_cols if c in adata_annot.obs.columns]
    annot_obs   = adata_annot.obs[available].copy()
    adata_annot.file.close()

    for col in available:
        adata.obs[col] = annot_obs.reindex(adata.obs_names)[col].values
    log.info("  Transferred: %s", available)

    # ── Step 7: UMAP plots ────────────────────────────────────────────────────
    color_cols = [c for c in [BATCH_KEY, "technology", "condition",
                               "leiden_scVI", "scanvi_Lake_label"]
                  if c in adata.obs.columns]
    if color_cols:
        fig = sc.pl.embedding(adata, basis="X_umap_scVI", color=color_cols,
                              ncols=2, show=False, return_fig=True)
        fig.savefig(os.path.join(PLOT_DIR, "umap_scvi_fullgene.pdf"),
                    bbox_inches="tight", dpi=150)
        plt.close(fig)

    # ── Step 8: Save ──────────────────────────────────────────────────────────
    log.info("Saving → %s", OUT_PATH)
    adata.write_h5ad(OUT_PATH)

    log.info("\n===== Complete =====")
    log.info("  Cells             : %d", adata.n_obs)
    log.info("  Genes (total)     : %d", adata.n_vars)
    log.info("  HVGs (model)      : %d", len(hvg_present))
    log.info("  scVI clusters     : %d", adata.obs["leiden_scVI"].nunique())
    if "leiden_scANVI" in adata.obs.columns:
        log.info("  scANVI clusters   : %d", adata.obs["leiden_scANVI"].nunique())
    log.info("  obsm keys         : %s", list(adata.obsm.keys()))
    log.info("  Output            : %s", OUT_PATH)


if __name__ == "__main__":
    main()
