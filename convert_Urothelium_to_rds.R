#!/usr/bin/env Rscript
# Convert Urothelium_cells.h5ad → Urothelium_cells.rds (Seurat object)
#
# Author: Xin Wang
# Date:   2026-04-23

suppressPackageStartupMessages({
  library(zellkonverter)
  library(Seurat)
  library(SingleCellExperiment)
})

base_dir <- "/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
h5ad_path <- file.path(base_dir, "integration_output", "Urothelium_cells.h5ad")
rds_path  <- file.path(base_dir, "integration_output", "Urothelium_cells.rds")

cat("Reading", h5ad_path, "\n")
sce <- readH5AD(h5ad_path, use_hdf5 = FALSE)
cat("  SCE dims:", nrow(sce), "genes x", ncol(sce), "cells\n")

cat("Converting SCE → Seurat ...\n")
so <- as.Seurat(sce, counts = "X", data = NULL)

# Transfer dimensional reductions if present
for (red in c("X_scVI", "X_scANVI", "X_umap", "X_umap_scVI", "X_umap_scANVI")) {
  if (red %in% reducedDimNames(sce)) {
    key <- sub("^X_", "", red)
    so[[key]] <- CreateDimReducObject(
      embeddings = reducedDim(sce, red),
      key        = paste0(key, "_"),
      assay      = DefaultAssay(so)
    )
    cat("  Added reduction:", key, "\n")
  }
}

cat("Saving to", rds_path, "\n")
saveRDS(so, rds_path)
cat("Done.\n")
