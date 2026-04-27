#!/usr/bin/env python3
"""
11_scvi_integrate_fullgene.py — scVI integration on the full-transcriptome
object (~25k genes, all 1.1M cells).

Key design
----------
  • Reads RenalUrothelium_allcells_fullgene.h5ad (built by 10_build_fullgene_h5ad.py).
  • Selects batch-aware HVGs (N_HVG, default 5 000) and marks them in var.
  • Trains scVI on the HVG subset — keeping all genes in the stored object.
  • Computes UMAP and Leiden clusters from the scVI latent space.
  • Transfers scanvi_Lake_label / scanvi_MKA_label annotations from the
    existing integrated annotated object (avoids re-running annotation).
  • Saves RenalUrothelium_allcells_scvi.h5ad with:
      layers["counts"]   — raw integer counts, all ~25k genes
      layers["lognorm"]  — log-normalised, all ~25k genes
      X                  — same as lognorm
      obsm["X_scVI"]     — 20-dim latent (batch-corrected)
      obsm["X_umap"]     — UMAP from scVI latent
      obs["leiden_scVI"] — Leiden clusters

Memory requirements
-------------------
  Loading the full-gene h5ad (~1.1M × 25k) needs ~200–400 G RAM.
  Run on a himem node with at least 500 G (1 300 G recommended).
  A GPU is strongly recommended for scVI training.

Usage
-----
  python 11_scvi_integrate_fullgene.py
  python 11_scvi_integrate_fullgene.py --n_hvg 8000 --n_epochs 200
"""

import argparse
import logging
import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR     = "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
SCRIPT_DIR   = os.path.join(BASE_DIR, "RenalUrotheliumScripts")
OUT_DIR      = os.path.join(SCRIPT_DIR, "output")
PLOT_DIR     = os.path.join(OUT_DIR, "plots")

IN_PATH      = os.path.join(OUT_DIR, "RenalUrothelium_allcells_fullgene.h5ad")
OUT_PATH     = os.path.join(OUT_DIR, "RenalUrothelium_allcells_scvi.h5ad")
MODEL_DIR    = os.path.join(OUT_DIR, "scvi_fullgene_model")
ANNOT_H5AD   = os.path.join(OUT_DIR, "RenalUrothelium_integrated_annotated.h5ad")

os.makedirs(PLOT_DIR, exist_ok=True)
os.makedirs(MODEL_DIR, exist_ok=True)

# ── Parameters ─────────────────────────────────────────────────────────────────
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

# ── Logging ────────────────────────────────────────────────────────────────────
logging.basicConfig(
    level  = logging.INFO,
    format = "%(asctime)s  %(levelname)s  %(message)s",
    datefmt= "%H:%M:%S",
)
log = logging.getLogger(__name__)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--n_hvg", type=int, default=5000,
                   help="Number of batch-aware HVGs for scVI training (default: 5000). "
                        "All other genes are retained in the object but not used for training.")
    p.add_argument("--n_epochs", type=int, default=N_EPOCHS,
                   help=f"Max training epochs (default: {N_EPOCHS})")
    p.add_argument("--skip_annotation_transfer", action="store_true",
                   help="Skip transferring Lake/MKA annotations from existing object")
    return p.parse_args()


def main():
    args   = parse_args()
    n_hvg  = args.n_hvg

    device = "cuda" if torch.cuda.is_available() else "cpu"
    log.info("PyTorch device: %s", device)
    if device == "cuda":
        log.info("  GPU: %s", torch.cuda.get_device_name(0))
    scvi.settings.dl_num_workers = int(os.environ.get("SLURM_CPUS_PER_TASK", 4))

    # ── Load full-gene object ─────────────────────────────────────────────────
    log.info("Loading %s ...", IN_PATH)
    adata = sc.read_h5ad(IN_PATH)
    log.info("  %d cells × %d genes", adata.n_obs, adata.n_vars)
    log.info("  Layers: %s", list(adata.layers.keys()))

    # ── Batch-aware HVG selection (mark only, do NOT subset) ─────────────────
    log.info("Selecting %d batch-aware HVGs (seurat_v3, batch=%s) ...", n_hvg, BATCH_KEY)
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes = n_hvg,
        batch_key   = BATCH_KEY,
        flavor      = "seurat_v3",
        layer       = "counts",
        subset      = False,          # ← mark HVGs but keep ALL genes
    )
    n_marked = int(adata.var["highly_variable"].sum())
    log.info("  HVGs marked: %d / %d total genes", n_marked, adata.n_vars)

    # ── Set up scVI on HVG subset ─────────────────────────────────────────────
    # Subset to HVGs for training — full adata keeps all genes
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    log.info("Training subset: %d cells × %d HVGs", adata_hvg.n_obs, adata_hvg.n_vars)

    cat_covs = [c for c in CAT_COVS if c in adata_hvg.obs.columns
                and adata_hvg.obs[c].nunique() > 1]

    scvi.model.SCVI.setup_anndata(
        adata_hvg,
        layer                      = "counts",
        batch_key                  = BATCH_KEY,
        categorical_covariate_keys = cat_covs or None,
    )

    model = scvi.model.SCVI(
        adata_hvg,
        n_latent        = N_LATENT,
        n_layers        = N_LAYERS,
        n_hidden        = N_HIDDEN,
        gene_likelihood = "nb",
        dispersion      = "gene-batch",
    )
    log.info("Model:\n%s", model)

    log.info("Training scVI (max_epochs=%d) ...", args.n_epochs)
    model.train(
        max_epochs              = args.n_epochs,
        batch_size              = BATCH_SIZE,
        early_stopping          = True,
        early_stopping_patience = 20,
        plan_kwargs             = {"lr": LR},
        accelerator             = "gpu" if device == "cuda" else "cpu",
    )

    # Save ELBO plot
    hist = model.history
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(hist["elbo_train"]["elbo_train"], label="train")
    if "elbo_validation" in hist:
        ax.plot(hist["elbo_validation"]["elbo_validation"], label="val")
    ax.set_xlabel("Epoch"); ax.set_ylabel("ELBO"); ax.legend()
    ax.set_title(f"scVI ELBO — full-gene ({n_hvg} HVGs)")
    fig.savefig(os.path.join(PLOT_DIR, "scvi_fullgene_elbo.pdf"), bbox_inches="tight")
    plt.close(fig)

    model.save(MODEL_DIR, overwrite=True)
    log.info("Model saved → %s", MODEL_DIR)

    # ── Latent representation (stored in full adata) ───────────────────────────
    log.info("Extracting scVI latent embedding ...")
    # get_latent_representation works on adata_hvg (same cells, HVG counts)
    adata.obsm["X_scVI"] = model.get_latent_representation()

    # ── Neighbours, UMAP, Leiden ──────────────────────────────────────────────
    log.info("Building neighbourhood graph ...")
    sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=N_NEIGHBOURS,
                    key_added="scVI")

    log.info("Leiden clustering (resolution=%.2f) ...", LEIDEN_RES)
    sc.tl.leiden(adata, neighbors_key="scVI", resolution=LEIDEN_RES,
                 key_added="leiden_scVI")

    log.info("Computing UMAP ...")
    sc.tl.umap(adata, neighbors_key="scVI")
    adata.obsm["X_umap_scVI"] = adata.obsm["X_umap"].copy()

    n_clusters = adata.obs["leiden_scVI"].nunique()
    log.info("  Clusters: %d", n_clusters)

    # ── Transfer existing annotations from integrated annotated object ─────────
    if not args.skip_annotation_transfer:
        annot_cols = ["scanvi_Lake_label", "scanvi_Lake_confidence",
                      "scanvi_MKA_label",  "scanvi_MKA_confidence"]
        log.info("Transferring annotations from %s ...", ANNOT_H5AD)
        adata_annot = sc.read_h5ad(ANNOT_H5AD, backed="r")
        annot_obs   = adata_annot.obs[
            [c for c in annot_cols if c in adata_annot.obs.columns]
        ].copy()
        adata_annot.file.close()

        shared = adata.obs_names.intersection(annot_obs.index)
        log.info("  Shared barcodes: %d / %d", len(shared), adata.n_obs)
        for col in annot_obs.columns:
            adata.obs[col] = annot_obs.loc[adata.obs_names, col] \
                             if len(shared) == adata.n_obs \
                             else annot_obs.reindex(adata.obs_names)[col]
        log.info("  Transferred: %s", list(annot_obs.columns))

    # ── UMAP plots ─────────────────────────────────────────────────────────────
    color_cols = [c for c in [BATCH_KEY, "technology", "condition",
                               "leiden_scVI", "scanvi_Lake_label"]
                  if c in adata.obs.columns]
    if color_cols:
        fig = sc.pl.embedding(adata, basis="X_umap_scVI", color=color_cols,
                              ncols=2, show=False, return_fig=True)
        fig.savefig(os.path.join(PLOT_DIR, "umap_scvi_fullgene.pdf"),
                    bbox_inches="tight", dpi=150)
        plt.close(fig)

    # ── Save ──────────────────────────────────────────────────────────────────
    log.info("Saving → %s", OUT_PATH)
    adata.write_h5ad(OUT_PATH)

    log.info("\n===== scVI full-gene integration complete =====")
    log.info("  Cells          : %d", adata.n_obs)
    log.info("  Genes (total)  : %d", adata.n_vars)
    log.info("  HVGs (trained) : %d", n_marked)
    log.info("  Clusters       : %d", n_clusters)
    log.info("  Output         : %s", OUT_PATH)


if __name__ == "__main__":
    main()
