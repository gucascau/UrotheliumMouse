################################################################################
# Script 04: Harmony Integration
#
# Resumes from the PCA checkpoint written by a prior run (Integrated_ScaledPCA.rds)
# and continues with:
#   5. RunHarmony  — batch correction by sample_id + technology
#   6. FindNeighbors / FindClusters
#   7. RunUMAP
#   8. Visualize
#   9. Save merged_harmony_integrated.rds
#
# Requires:
#   integration_output/Integrated_ScaledPCA.rds   (ScaleData + RunPCA already done)
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

OUT_PATH  <- file.path(OUT_DIR, "merged_harmony_integrated.rds")

# ── Parameters ────────────────────────────────────────────────────────────────
N_PCS        <- 50
HARMONY_VARS <- c("sample_id", "technology")
HARMONY_DIMS <- 1:30
RESOLUTION   <- 0.5


################################################################################
# STEP 1-4: Load PCA checkpoint (ScaleData + PCA already done)
################################################################################

PCA_PATH <- file.path(OUT_DIR, "Integrated_ScaledPCA.rds")
if (!file.exists(PCA_PATH))
  stop("Checkpoint not found: ", PCA_PATH,
       "\nRun steps 1-4 first (ScaleData + RunPCA).")

message("Loading PCA checkpoint...")
merged <- readRDS(PCA_PATH)
message(sprintf("  Loaded: %d cells × %d genes", ncol(merged), nrow(merged)))
message(sprintf("  PCA dims: %s", paste(dim(Embeddings(merged, "pca")), collapse = " × ")))
log_mem("after checkpoint load")

# Drop scale.data if still present to free RAM before Harmony
if (!is.null(merged[["RNA"]]$scale.data)) {
  merged[["RNA"]][["scale.data"]] <- NULL
  gc()
  log_mem("after scale.data freed")
}


################################################################################
# STEP 5: RunHarmony  (skipped if checkpoint exists)
################################################################################

HARMONY_CHECKPOINT <- file.path(OUT_DIR, "harmony_embeddings.rds")

if (file.exists(HARMONY_CHECKPOINT)) {
  message("Harmony checkpoint found — loading embeddings, skipping RunHarmony...")
  harmony_mat <- readRDS(HARMONY_CHECKPOINT)
  merged[["harmony"]] <- CreateDimReducObject(
    embeddings = harmony_mat,
    key        = "harmony_",
    assay      = DefaultAssay(merged)
  )
  rm(harmony_mat); gc()
  log_mem("after loading Harmony checkpoint")
} else {
  message(sprintf("RunHarmony (batch = %s)...",
                  paste(HARMONY_VARS, collapse = " + ")))
  pca_mat <- Embeddings(merged, "pca")
  harmony_mat <- suppressWarnings(
    harmony::HarmonyMatrix(
      data_mat          = pca_mat,
      meta_data         = merged@meta.data,
      vars_use          = HARMONY_VARS,
      do_pca            = FALSE,
      max.iter.harmony  = 20,
      verbose           = FALSE
    )
  )
  merged[["harmony"]] <- CreateDimReducObject(
    embeddings = harmony_mat,
    key        = "harmony_",
    assay      = DefaultAssay(merged)
  )
  rm(pca_mat, harmony_mat); gc()
  message("  Harmony done")
  log_mem("after Harmony")

  saveRDS(Embeddings(merged, "harmony"), HARMONY_CHECKPOINT)
  message("  Harmony embedding checkpoint saved")
}

# Replace the RNA assay with a 1-gene placeholder to free all expression data
# (~80-100 GB) before building the annoy index for 2.8M cells.
# Seurat v5 Assay5 forbids removing its last layer, so rebuilding is the only
# way to guarantee a full drop of counts + data + scale.data in one step.
cell_names  <- colnames(merged)
empty_mat   <- Matrix::sparseMatrix(
  i = 1L, j = seq_along(cell_names),
  x = 0, dims = c(1L, length(cell_names)),
  dimnames = list("placeholder", cell_names)
)
merged[["RNA"]] <- CreateAssay5Object(counts = empty_mat)
rm(empty_mat, cell_names)
merged[["pca"]] <- NULL
gc()
log_mem("after stripping assay data for FindNeighbors")

################################################################################
# STEP 6: FindNeighbors + FindClusters
################################################################################

message("FindNeighbors (annoy, k=5)...")
merged <- FindNeighbors(
  merged,
  reduction = "harmony",
  dims      = HARMONY_DIMS,
  nn.method = "annoy",
  k.param   = 5,
  annoy.metric = "euclidean",
  n.trees   = 50,
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
