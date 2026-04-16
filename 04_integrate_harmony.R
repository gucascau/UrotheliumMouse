################################################################################
# Script 04: Harmony Integration
#
# Reads the merged + log-normalized object written by 03_integrate_harmony.R
# (merged_normalized.rds) and the companion HVG list (hvg_list.rds), then runs:
#   1. Restore HVG selection
#   2. JoinLayers  — collapse per-sample layers into one matrix
#   3. ScaleData   — regress out pct_mt
#   4. RunPCA      — 50 PCs on HVGs
#   5. RunHarmony  — batch correction by sample_id + technology
#   6. FindNeighbors / FindClusters
#   7. RunUMAP
#   8. Visualize
#   9. Save merged_harmony_integrated.rds
#
# Run AFTER 03_integrate_harmony.R has written:
#   integration_output/merged_normalized.rds
#   integration_output/hvg_list.rds
################################################################################

library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)

# FindNeighbors uses future for parallelism; on 2.8M cells the SNN index and
# query matrix exceed the 500 MiB default.  Raise the limit to 8 GiB.
options(future.globals.maxSize = 8 * 1024^3)

# Helper: log current R memory usage at key steps
log_mem <- function(label) {
  g  <- gc(verbose = FALSE)
  mb <- sum(g[, 2])
  message(sprintf("  [mem] %s: %.1f GB in use", label, mb / 1024))
}

# ── Paths ─────────────────────────────────────────────────────────────────────
DATA_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
OUT_DIR  <- file.path(DATA_DIR, "integration_output")
dir.create(OUT_DIR, showWarnings = FALSE)

NORM_PATH <- file.path(OUT_DIR, "merged_normalized.rds")
HVG_PATH  <- file.path(OUT_DIR, "hvg_list.rds")
OUT_PATH  <- file.path(OUT_DIR, "merged_harmony_integrated.rds")

# ── Parameters ────────────────────────────────────────────────────────────────
N_PCS        <- 50
HARMONY_VARS <- c("sample_id", "technology")
HARMONY_DIMS <- 1:30
RESOLUTION   <- 0.5


################################################################################
# STEP 1: Load merged normalized object + HVG list
################################################################################

for (p in c(NORM_PATH, HVG_PATH)) {
  if (!file.exists(p))
    stop("Required file not found: ", p,
         "\nRun 03_integrate_harmony.R first.")
}

message("Loading merged normalized object...")
merged <- readRDS(NORM_PATH)
message(sprintf("  Loaded: %d cells × %d genes", ncol(merged), nrow(merged)))
message(sprintf("  Layers: %s", paste(Layers(merged), collapse = ", ")))
log_mem("after load")

message("Loading HVG list...")
hvg_obj <- readRDS(HVG_PATH)
VariableFeatures(merged) <- hvg_obj$hvg
message(sprintf("  HVGs restored: %d", length(VariableFeatures(merged))))
rm(hvg_obj); gc()


################################################################################
# STEP 2: JoinLayers
################################################################################

message("JoinLayers — collapsing per-sample layers...")
merged <- JoinLayers(merged)
message(sprintf("  Layers after join: %s",
                paste(Layers(merged), collapse = ", ")))
log_mem("after JoinLayers")


################################################################################
# STEP 3: ScaleData
################################################################################

message("ScaleData (regress pct_mt, HVGs only)...")
merged <- ScaleData(merged,
                    features        = VariableFeatures(merged),
                    vars.to.regress = "pct_mt",
                    verbose         = FALSE)
message("  ScaleData done")
log_mem("after ScaleData")


################################################################################
# STEP 4: RunPCA
################################################################################

message(sprintf("RunPCA (%d PCs)...", N_PCS))
merged <- RunPCA(merged, npcs = N_PCS, verbose = FALSE)
message("  PCA done")

# save into the Seurat object for later checkpointing — allows re-running FindNeighbors/UMAP
saveRDS(merged,
        file.path(OUT_DIR, "Integrated_ScaledPCA.rds"))

# scale.data is dense (HVGs × 2.8M cells, ~69 GB) — only needed for RunPCA.
# Drop it immediately to free RAM before Harmony + FindNeighbors.
merged[["RNA"]]$scale.data <- NULL
gc()
log_mem("after PCA + scale.data freed")


################################################################################
# STEP 5: RunHarmony
################################################################################

message(sprintf("RunHarmony (batch = %s)...",
                paste(HARMONY_VARS, collapse = " + ")))
# suppressWarnings: silences k-means "Quick-TRANSfer" non-convergence warnings
# at 2.8M cells — Harmony still converges correctly.
suppressWarnings(
  merged <- RunHarmony(
    object         = merged,
    group.by.vars  = HARMONY_VARS,
    reduction      = "pca",
    reduction.save = "harmony",
    max_iter       = 20,
    verbose        = FALSE
  )
)
message("  Harmony done")
log_mem("after Harmony")

# Checkpoint: save Harmony embedding so FindNeighbors/UMAP can be rerun
# without repeating the 4-hour ScaleData + PCA + Harmony steps.
saveRDS(Embeddings(merged, "harmony"),
        file.path(OUT_DIR, "harmony_embeddings.rds"))
message("  Harmony embedding checkpoint saved")

################################################################################
# STEP 6: FindNeighbors + FindClusters
################################################################################

message("FindNeighbors (annoy, k=15)...")
# annoy: approximate NN — far lower peak memory than exact NN for 2.8M cells
merged <- FindNeighbors(
  merged,
  reduction = "harmony",
  dims      = HARMONY_DIMS,
  nn.method = "annoy",
  k.param   = 15,
  verbose   = FALSE
)
log_mem("after FindNeighbors")

message(sprintf("FindClusters (resolution = %.2f)...", RESOLUTION))
merged <- FindClusters(merged, resolution = RESOLUTION, verbose = FALSE)
message(sprintf("  Clusters: %d", length(unique(merged$seurat_clusters))))
gc()


################################################################################
# STEP 7: RunUMAP
################################################################################

message("RunUMAP on Harmony embedding...")
merged <- RunUMAP(merged, reduction = "harmony", dims = HARMONY_DIMS,
                  reduction.name = "umap", verbose = FALSE)
message("  UMAP done")


################################################################################
# STEP 8: Visualize
################################################################################

message("Generating UMAP plots...")

p1 <- DimPlot(merged, group.by = "sample_id",      reduction = "umap",
              label = FALSE, raster = TRUE) + ggtitle("By Sample (batch)")
p2 <- DimPlot(merged, group.by = "condition",       reduction = "umap",
              label = FALSE, raster = TRUE) + ggtitle("By Condition")
p3 <- DimPlot(merged, group.by = "technology",      reduction = "umap",
              label = FALSE, raster = TRUE) + ggtitle("By Technology")
p4 <- DimPlot(merged, group.by = "seurat_clusters", reduction = "umap",
              label = TRUE,  raster = TRUE) + ggtitle("Clusters")
p5 <- DimPlot(merged, group.by = "paper",
              reduction = "umap",
              label = FALSE, raster = TRUE) + ggtitle("By Paper")
p6 <- DimPlot(merged, group.by = "gsm_id",          reduction = "umap",
              label = FALSE, raster = TRUE) + ggtitle("By GSM/GSE ID")

pdf(file.path(OUT_DIR, "UMAP_overview.pdf"), width = 18, height = 24)
print(p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 2))
dev.off()

if ("cell_type_original" %in% colnames(merged@meta.data)) {
  p7 <- DimPlot(merged, group.by = "cell_type_original",
                reduction = "umap",
                label = TRUE, repel = TRUE, raster = TRUE) +
        ggtitle("Known cell type annotations")
  pdf(file.path(OUT_DIR, "UMAP_known_annotations.pdf"), width = 12, height = 10)
  print(p7)
  dev.off()
}

# Harmony convergence plot (if available)
harmony_obj <- merged[["harmony"]]
conv_plot   <- tryCatch(harmony_obj@misc$convergence_plot, error = function(e) NULL)
if (!is.null(conv_plot)) {
  pdf(file.path(OUT_DIR, "harmony_convergence.pdf"), width = 6, height = 4)
  print(conv_plot)
  dev.off()
}


################################################################################
# STEP 9: Save
################################################################################

message(sprintf("Saving integrated object to %s ...", OUT_PATH))
saveRDS(merged, OUT_PATH)

message("\n===== Harmony integration complete =====")
message(sprintf("Total cells: %s", format(ncol(merged), big.mark = ",")))
message(sprintf("Clusters:    %d", length(unique(merged$seurat_clusters))))
message(sprintf("Batches:     %d samples", length(unique(merged$sample_id))))
message(sprintf("Output:      %s", OUT_PATH))
