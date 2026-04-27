#!/usr/bin/env Rscript
# 07b_convert_annotated_to_rds.R — Convert RenalUrothelium_integrated_annotated.h5ad
# to a Seurat RDS object, transferring all reductions and annotation columns.
#
# Input  : output/RenalUrothelium_integrated_annotated.h5ad
# Output : output/RenalUrothelium_integrated_annotated.rds

suppressPackageStartupMessages({
  library(zellkonverter)
  library(Seurat)
  library(SingleCellExperiment)
})

base_dir  <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
out_dir   <- file.path(base_dir, "RenalUrotheliumScripts", "output")
h5ad_path <- file.path(out_dir, "RenalUrothelium_integrated_annotated.h5ad")
rds_path  <- file.path(out_dir, "RenalUrothelium_integrated_annotated.rds")

cat("Reading", h5ad_path, "\n")
sce <- readH5AD(h5ad_path, use_hdf5 = FALSE)
cat("  SCE dims:", nrow(sce), "genes ×", ncol(sce), "cells\n")

cat("Converting SCE → Seurat ...\n")
# counts = raw integer counts (layers["counts"] in h5ad)
# data   = log-normalised expression (X in h5ad, i.e. normalize_total + log1p)
# This allows pseudobulk DEG (DESeq2/edgeR) via counts and FindMarkers via data.
so <- as.Seurat(sce, counts = "counts", data = "X")

# Transfer dimensional reductions
for (red in c("X_scVI", "X_scANVI", "X_umap_scVI", "X_umap_scANVI", "X_umap")) {
  if (red %in% reducedDimNames(sce)) {
    key    <- sub("^X_", "", red)
    emb    <- reducedDim(sce, red)
    so[[key]] <- CreateDimReducObject(
      embeddings = emb,
      key        = paste0(gsub("_", "", key), "_"),
      assay      = DefaultAssay(so)
    )
    cat("  Added reduction:", key, "\n")
  }
}

# Report annotation columns
for (col in c("scanvi_Lake_label", "scanvi_Lake_confidence",
              "scanvi_MKA_label",  "scanvi_MKA_confidence")) {
  if (col %in% colnames(so@meta.data)) {
    cat(sprintf("  %s: %d non-NA values\n", col,
                sum(!is.na(so@meta.data[[col]]))))
  }
}

cat("Saving to", rds_path, "\n")
saveRDS(so, rds_path)
cat("Done.\n")
cat("  Cells:", ncol(so), "\n")
cat("  Genes:", nrow(so), "\n")
