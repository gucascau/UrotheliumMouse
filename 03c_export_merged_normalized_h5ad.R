################################################################################
# Script 03c: Export merged_normalized.rds → AnnData (.h5ad)
#
# Reads the merged + log-normalised Seurat v5 object written by
# 03_integrate_harmony.R (merged_normalized.rds) and exports it as an h5ad
# file for downstream Python / scVI / scanpy workflows.
#
# AnnData layout:
#   adata.X              — log-normalised data (sparse float, cells × genes)
#   adata.layers["counts"] — raw integer counts, if present in merged object
#   adata.obs            — cell metadata (sample_id, condition, technology, …)
#   adata.var            — gene metadata + highly_variable flag (from hvg_list.rds)
#   adata.uns["hvg"]     — list of HVG gene names
#
# Memory: loading merged_normalized.rds + JoinLayers + transpose requires
#         ~400–512 GB RAM.  Run on a himem node (512 GB).
#
# Run AFTER 03_integrate_harmony.R has written:
#   integration_output/merged_normalized.rds
#   integration_output/hvg_list.rds
################################################################################

library(Seurat)
library(Matrix)
library(reticulate)
library(zellkonverter)

# ── Paths ─────────────────────────────────────────────────────────────────────
DATA_DIR  <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
OUT_DIR   <- file.path(DATA_DIR, "integration_output")
dir.create(OUT_DIR, showWarnings = FALSE)

NORM_PATH <- file.path(OUT_DIR, "merged_normalized.rds")
HVG_PATH  <- file.path(OUT_DIR, "hvg_list.rds")
H5AD_PATH <- file.path(OUT_DIR, "merged_normalized.h5ad")

message("=== 03c: Export merged_normalized.rds → h5ad ===")
message(sprintf("Input  : %s", NORM_PATH))
message(sprintf("Output : %s", H5AD_PATH))

for (p in c(NORM_PATH, HVG_PATH)) {
  if (!file.exists(p))
    stop(sprintf("Required file not found: %s\nRun 03_integrate_harmony.R first.", p))
}


################################################################################
# STEP 1: Load HVG list
################################################################################

message("\n[1/5] Loading HVG list...")
hvg_data  <- readRDS(HVG_PATH)
hvg_genes <- hvg_data$hvg
all_genes  <- hvg_data$all_genes
message(sprintf("  HVGs: %d / total genes: %d", length(hvg_genes), length(all_genes)))


################################################################################
# STEP 2: Load merged normalised Seurat object
################################################################################

message("\n[2/5] Loading merged_normalized.rds (this may take several minutes)...")
merged <- readRDS(NORM_PATH)
message(sprintf("  Loaded: %d cells × %d genes", ncol(merged), nrow(merged)))
message(sprintf("  Layers present: %s", paste(Layers(merged), collapse = ", ")))


################################################################################
# STEP 3: JoinLayers — collapse per-sample layers into single matrices
################################################################################

message("\n[3/5] JoinLayers — collapsing per-sample layers...")
DefaultAssay(merged) <- "RNA"
merged <- JoinLayers(merged)
message(sprintf("  Layers after join: %s", paste(Layers(merged), collapse = ", ")))


################################################################################
# STEP 4: Export to h5ad for scVI integration or other method integration
################################################################################

# Keep only data + counts before converting.  Setting a layer to NULL in
# Seurat v5 leaves an empty stub that trips SCE dimension validation
# ("all assays must have the same nrow and ncol").  DietSeurat() removes
# layers cleanly without leaving stubs.
message("\n[4/5] Trimming to data + counts layers (removing scale.data stub)...")
merged <- DietSeurat(merged, layers = c("data", "counts"), assays = "RNA")
message(sprintf("  Layers for export: %s", paste(Layers(merged), collapse = ", ")))

message("\n[4/5] Converting to SingleCellExperiment...")
sce <- as.SingleCellExperiment(merged)

message("\n[5/5] Writing h5ad file...")
writeH5AD(
  sce,
  file = H5AD_PATH
)
cat("h5ad written to:", H5AD_PATH, "\n")
