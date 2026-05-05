#!/usr/bin/env Rscript
# 04b_convert_renal_epithelial_to_rds.R — Convert
# RenalUrothelium_renal_epithelial.h5ad → Seurat RDS.
#
# The h5ad stores:
#   X              : log-normalised counts (float32)
#   layers/counts  : raw integer counts
#   obs            : includes epithelial_compartment, assign_method

suppressPackageStartupMessages({
  library(zellkonverter)
  library(Seurat)
  library(SingleCellExperiment)
})

base_dir  <- paste0("/vast0/home/gdjacksonlab/lab/xxw004/UUO/",
                    "Datasets/Mouse/UsedSingleCells")
out_dir   <- file.path(base_dir, "RenalUrotheliumScripts", "output")
h5ad_path <- file.path(out_dir, "RenalUrothelium_renal_epithelial.h5ad")
rds_path  <- file.path(out_dir, "RenalUrothelium_renal_epithelial.rds")

cat("Reading", h5ad_path, "\n")
sce <- readH5AD(h5ad_path, use_hdf5 = FALSE)
cat("  SCE dims:", nrow(sce), "genes ×", ncol(sce), "cells\n")

# X is log-norm; use it as the Seurat 'counts' slot (same convention as
# 04_convert_to_rds.R, so 05b_harmony_recluster can treat it identically).
cat("Converting SCE → Seurat ...\n")
so <- as.Seurat(sce, counts = "X", data = NULL)

# as.Seurat may name the assay "originalexp" or "X" rather than "RNA";
# rename it so downstream scripts can rely on so[["RNA"]].
default_assay <- DefaultAssay(so)
if (default_assay != "RNA") {
  so <- RenameAssays(so, assay = default_assay, new.assay.name = "RNA")
  cat(sprintf("  Renamed assay '%s' → 'RNA'\n", default_assay))
}

# Transfer any dimensional reductions (may not exist pre-integration)
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

# Report compartment composition
if ("epithelial_compartment" %in% colnames(so@meta.data)) {
  cat("\nCells per epithelial_compartment:\n")
  print(sort(table(so@meta.data$epithelial_compartment), decreasing = TRUE))
}
if ("assign_method" %in% colnames(so@meta.data)) {
  cat("\nCells per assign_method:\n")
  print(table(so@meta.data$assign_method))
}

cat("\nSaving to", rds_path, "\n")
saveRDS(so, rds_path)
cat("Done.\n")
cat("  Cells:", ncol(so), "\n")
cat("  Genes:", nrow(so), "\n")
