#!/usr/bin/env Rscript
# 10b_convert_fullgene_to_rds.R — Convert the full-gene h5ad (built by
# 10_build_fullgene_h5ad.py) to a Seurat RDS, then transfer the scVI/scANVI
# UMAP reductions from the integrated annotated RDS.
#
# Input
# -----
#   output/RenalUrothelium_allcells_fullgene.h5ad  (1.1M cells × ~25k genes)
#   output/RenalUrothelium_integrated_annotated.rds    (for UMAP reductions)
#
# Output
# ------
#   output/RenalUrothelium_allcells_fullgene.rds
#     counts slot  — raw integer counts (~25k genes)
#     data   slot  — log-normalised expression
#     All metadata columns from integrated object
#     scVI, scANVI, umap_scVI, umap_scANVI reductions

suppressPackageStartupMessages({
  library(zellkonverter)
  library(Seurat)
  library(SingleCellExperiment)
})

BASE_DIR   <- paste0("/vast0/home/gdjacksonlab/lab/xxw004/UUO/",
                     "Datasets/Mouse/UsedSingleCells")
SCRIPT_DIR <- file.path(BASE_DIR, "RenalUrotheliumScripts")
OUT_DIR    <- file.path(SCRIPT_DIR, "output")

H5AD_IN    <- file.path(OUT_DIR, "RenalUrothelium_allcells_fullgene.h5ad")
INT_RDS    <- file.path(OUT_DIR, "RenalUrothelium_integrated_annotated.rds")
RDS_OUT    <- file.path(OUT_DIR, "RenalUrothelium_allcells_fullgene.rds")

# ── Read full-gene h5ad ───────────────────────────────────────────────────────
cat("Reading", H5AD_IN, "\n")
sce <- readH5AD(H5AD_IN, use_hdf5 = FALSE)
cat("  SCE dims:", nrow(sce), "genes ×", ncol(sce), "cells\n")
cat("  Assays:", assayNames(sce), "\n")

# ── Convert to Seurat ─────────────────────────────────────────────────────────
cat("Converting SCE → Seurat ...\n")
so <- as.Seurat(sce, counts = "counts", data = "lognorm")
rm(sce)
gc()
cat("  Genes:", nrow(so), "\n")
cat("  Cells:", ncol(so), "\n")

# ── Transfer reductions from integrated RDS ───────────────────────────────────
cat("\nLoading integrated RDS for reductions:", INT_RDS, "\n")
so_int <- readRDS(INT_RDS)

# Keep only cells present in full-gene object
shared_cells <- intersect(colnames(so), colnames(so_int))
cat("  Shared cells:", length(shared_cells), "\n")
so_int <- so_int[, shared_cells]
so     <- so[,    shared_cells]

for (red in c("scVI", "scANVI", "umap_scVI", "umap_scANVI", "umap")) {
  if (red %in% names(so_int@reductions)) {
    emb    <- Embeddings(so_int, red)
    so[[red]] <- CreateDimReducObject(
      embeddings = emb,
      key        = paste0(gsub("_", "", red), "_"),
      assay      = DefaultAssay(so)
    )
    cat("  Transferred reduction:", red, "\n")
  }
}
rm(so_int)
gc()

# ── Save ──────────────────────────────────────────────────────────────────────
cat("\nSaving →", RDS_OUT, "\n")
saveRDS(so, RDS_OUT)
cat("Done.\n")
cat("  Cells:", ncol(so), "\n")
cat("  Genes:", nrow(so), "\n")
cat("  Reductions:", names(so@reductions), "\n")
cat("  Meta columns:", ncol(so@meta.data), "\n")
