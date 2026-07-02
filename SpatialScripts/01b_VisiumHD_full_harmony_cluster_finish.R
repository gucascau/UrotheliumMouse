################################################################################
# 01b_VisiumHD_full_harmony_cluster_finish.R
#
# VisiumHD_harmony_integrated.rds (bladder + 2 kidneys, 1,264,551 cells) was
# saved after merge/normalize/scale/sketch/PCA-on-sketch, but the pipeline
# was interrupted before Harmony batch-correction, clustering, UMAP, and
# ProjectData — those steps were only completed for the kidney-only pivot.
# This script finishes the full 3-dataset object at that same point,
# mirroring 01_VisiumHD_integrate_harmony_sketch.R steps 7-10.
#
# Input : output/VisiumHD_harmony_integrated.rds  (dataset = Visium_HD_Bladder,
#           Visium_HD_3prime_Kidney, Visium_HD_Kidney; sketch = 5000 cells/ds;
#           pca.sketch already computed)
# Output: output/VisiumHD_harmony_integrated_clustered.rds
#         output/VisiumHD_harmony_full_clusters.pdf
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

options(future.globals.maxSize = 16 * 1024^3)

log_mem <- function(label) {
  g  <- gc(verbose = FALSE)
  mb <- sum(g[, 2])
  message(sprintf("  [mem] %s: %.1f GB in use", label, mb / 1024))
}

OUT_DIR  <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/SpatialScripts/output"
IN_RDS   <- file.path(OUT_DIR, "VisiumHD_harmony_integrated.rds")
OUT_RDS  <- file.path(OUT_DIR, "VisiumHD_harmony_integrated_clustered.rds")
OUT_PDF  <- file.path(OUT_DIR, "VisiumHD_harmony_full_clusters.pdf")

ASSAY <- "Spatial.008um"

message("==> Loading: ", IN_RDS)
object <- readRDS(IN_RDS)
message(sprintf("  %d cells x %d genes, datasets: %s",
  ncol(object), nrow(object), paste(unique(object$dataset), collapse = ", ")))
log_mem("after load")

DefaultAssay(object) <- "sketch"

# ── Step 6 (already done at load, verify) ──────────────────────────────────────
if (!"pca.sketch" %in% names(object@reductions)) {
  message("==> PCA on sketch (not found, computing) ...")
  object <- FindVariableFeatures(object, assay = "sketch", nfeatures = 3000)
  object <- ScaleData(object, assay = "sketch")
  object <- RunPCA(object, assay = "sketch", reduction.name = "pca.sketch", npcs = 30)
}

# ── Step 7: Harmony batch correction on sketch ─────────────────────────────────
message("==> Running Harmony on sketch (batch = dataset) ...")
object <- RunHarmony(
  object          = object,
  group.by.vars   = "dataset",
  reduction       = "pca.sketch",
  reduction.save  = "harmony.sketch",
  theta           = 2,
  max_iter        = 20,
  verbose         = TRUE
)
log_mem("after Harmony")

# ── Step 8: Cluster and embed sketched cells (harmony embedding) ──────────────
message("==> Clustering sketched data ...")
object <- FindNeighbors(object,
  reduction = "harmony.sketch",
  dims      = 1:30,
  assay     = "sketch"
)
object <- FindClusters(object,
  resolution   = 0.5,
  cluster.name = "seurat_cluster.harmony.sketched"
)
object <- RunUMAP(object,
  reduction      = "harmony.sketch",
  dims           = 1:30,
  reduction.name = "umap.harmony.sketch",
  return.model   = TRUE
)
log_mem("after sketch clustering")

# ── Step 9: Project back to full dataset ───────────────────────────────────────
message("==> Projecting clusters back to full dataset ...")
object <- ProjectData(
  object             = object,
  assay              = ASSAY,
  full.reduction     = "full.pca.sketch",
  sketched.assay     = "sketch",
  sketched.reduction = "harmony.sketch",
  umap.model         = "umap.harmony.sketch",
  dims               = 1:30,
  refdata            = list(
    seurat_cluster.harmony.projected = "seurat_cluster.harmony.sketched"
  )
)
log_mem("after projection")

message(sprintf("  %d clusters found", n_distinct(object$seurat_cluster.harmony.projected)))
print(table(object$seurat_cluster.harmony.projected, object$dataset))

# ── Step 10: Visualize and save ─────────────────────────────────────────────────
message("==> Saving plots ...")
Idents(object) <- "seurat_cluster.harmony.projected"

pdf(OUT_PDF, width = 14, height = 6)

p_umap <- DimPlot(object,
  reduction = "umap.harmony.sketch",
  group.by  = "seurat_cluster.harmony.projected",
  label     = TRUE,
  repel     = TRUE
) + ggtitle("Harmony+Sketch clusters (UMAP) — bladder + kidneys")

p_batch <- DimPlot(object,
  reduction = "umap.harmony.sketch",
  group.by  = "dataset",
  label     = FALSE
) + ggtitle("Dataset identity")

print(p_umap | p_batch)

for (samp in unique(object$dataset)) {
  cells_in_samp <- WhichCells(object, expression = dataset == samp)
  sub_obj <- subset(object, cells = cells_in_samp)
  Idents(sub_obj) <- "seurat_cluster.harmony.projected"
  tryCatch({
    print(
      SpatialDimPlot(sub_obj, label = FALSE, pt.size.factor = 3) +
        ggtitle(paste("Spatial clusters:", samp))
    )
  }, error = function(e) {
    message("  SpatialDimPlot failed for ", samp, ": ", e$message)
  })
}

dev.off()
message("  Saved: ", OUT_PDF)

message("==> Saving RDS ...")
saveRDS(object, OUT_RDS)
message("  Saved: ", OUT_RDS)

message("==> Done.")
log_mem("final")
