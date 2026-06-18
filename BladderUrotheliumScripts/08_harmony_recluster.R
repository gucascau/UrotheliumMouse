################################################################################
# 08_harmony_recluster.R — Harmony reclustering of bladder urothelial cells
#
# Input : output/BladderUrothelium_uro_cells_scvi.rds
#         counts layer = raw counts
#         data layer   = log-normalised expression (layers["lognorm"])
#
# Steps:
#   1. Load Seurat object
#   2. Verify counts/data layers
#   3. Compute pct_mt + FindVariableFeatures (3 000 HVGs)
#   4. ScaleData (regress pct_mt)
#   5. RunPCA (20 PCs)
#   6. RunHarmony (sample_id; +technology if >1 level)
#   7. FindNeighbors + FindClusters (res = 0.5)
#   8. RunUMAP
#   9. UMAP + marker plots → PDFs
#  10. Save output/BladderUrothelium_uro_cells_harmony_integrated.rds
################################################################################

library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)

options(future.globals.maxSize = 4 * 1024^3)

log_mem <- function(label) {
  g  <- gc(verbose = FALSE)
  mb <- sum(g[, 2])
  message(sprintf("  [mem] %s: %.1f GB in use", label, mb / 1024))
}

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
OUT_DIR  <- file.path(BASE_DIR, "BladderUrothelium", "BladderUrotheliumScripts", "output")
dir.create(OUT_DIR, showWarnings = FALSE)

IN_PATH  <- file.path(OUT_DIR, "BladderUrothelium_uro_cells_scvi.rds")
OUT_PATH <- file.path(OUT_DIR, "BladderUrothelium_uro_cells_harmony_integrated.rds")

# ── Parameters ────────────────────────────────────────────────────────────────
N_HVG        <- 3000
N_PCS        <- 20
HARMONY_DIMS <- 1:20
RESOLUTION   <- 0.5


################################################################################
# STEP 1: Load
################################################################################

if (!file.exists(IN_PATH))
  stop("Input not found: ", IN_PATH)

message("Loading BladderUrothelium_uro_cells_scvi.rds ...")
so <- readRDS(IN_PATH)
message(sprintf("  Loaded: %d cells × %d genes", ncol(so), nrow(so)))
log_mem("after load")

for (col in c("sample_id", "condition", "technology")) {
  if (col %in% colnames(so@meta.data)) {
    message(sprintf("  Cells per %s:", col))
    print(sort(table(so@meta.data[[col]]), decreasing = TRUE))
  }
}


################################################################################
# STEP 2: Verify expression layers
################################################################################

message("Verifying RNA assay layers ...")
rna_layers <- Layers(so[["RNA"]])
message("  RNA layers: ", paste(rna_layers, collapse = ", "))

for (layer in c("counts", "data")) {
  if (!layer %in% rna_layers)
    stop("RNA assay is missing required layer: ", layer)
}
log_mem("after layer check")


################################################################################
# STEP 3: pct_mt + FindVariableFeatures
################################################################################

message("Computing pct_mt ...")
so[["pct_mt"]] <- PercentageFeatureSet(so, pattern = "^mt-")
message(sprintf("  pct_mt range: %.2f – %.2f", min(so$pct_mt), max(so$pct_mt)))

n_na <- sum(is.na(so$pct_mt))
if (n_na > 0) {
  med_val <- median(so$pct_mt, na.rm = TRUE)
  so$pct_mt[is.na(so$pct_mt)] <- med_val
  message(sprintf("  Imputed %d NAs in pct_mt with median (%.4f)", n_na, med_val))
}

message(sprintf("FindVariableFeatures (nfeatures = %d) ...", N_HVG))
so <- FindVariableFeatures(so, selection.method = "vst",
                           nfeatures = N_HVG, verbose = FALSE)
message(sprintf("  HVGs: %d", length(VariableFeatures(so))))


################################################################################
# STEP 4: ScaleData
################################################################################

message("ScaleData (regress pct_mt, HVGs only) ...")
so <- ScaleData(so,
                features        = VariableFeatures(so),
                vars.to.regress = "pct_mt",
                verbose         = FALSE)
message("  ScaleData done")
log_mem("after ScaleData")


################################################################################
# STEP 5: RunPCA
################################################################################

message(sprintf("RunPCA (%d PCs) ...", N_PCS))
so <- RunPCA(so, npcs = N_PCS, verbose = FALSE)
message("  PCA done")
log_mem("after PCA")

so[["RNA"]]$counts     <- NULL
so[["RNA"]]$data       <- NULL
so[["RNA"]]$scale.data <- NULL
gc()
log_mem("after freeing expression data")


################################################################################
# STEP 6: RunHarmony
################################################################################

harmony_vars <- "sample_id"
for (v in c("technology", "source")) {
  if (v %in% colnames(so@meta.data) &&
      length(unique(so@meta.data[[v]])) > 1) {
    harmony_vars <- c(harmony_vars, v)
  }
}
message(sprintf("RunHarmony (batch = %s) ...", paste(harmony_vars, collapse = " + ")))

so <- RunHarmony(
  so,
  group.by.vars    = harmony_vars,
  reduction        = "pca",
  reduction.save   = "harmony",
  plot_convergence = FALSE,
  verbose          = FALSE
)
so[["pca"]] <- NULL
gc()
message("  Harmony done")
log_mem("after Harmony")


################################################################################
# STEP 7: FindNeighbors + FindClusters
################################################################################

message("FindNeighbors (annoy, k=20) ...")
so <- FindNeighbors(
  so,
  reduction    = "harmony",
  dims         = HARMONY_DIMS,
  nn.method    = "annoy",
  k.param      = 20,
  annoy.metric = "euclidean",
  n.trees      = 50,
  verbose      = FALSE
)
log_mem("after FindNeighbors")

message(sprintf("FindClusters (resolution = %.2f) ...", RESOLUTION))
so <- FindClusters(so, resolution = RESOLUTION, verbose = FALSE)
message(sprintf("  Clusters: %d", length(unique(so$seurat_clusters))))
gc()


################################################################################
# STEP 8: RunUMAP
################################################################################

message("RunUMAP on Harmony embedding ...")
so <- RunUMAP(so, reduction = "harmony", dims = HARMONY_DIMS,
              reduction.name = "umap_harmony", verbose = FALSE)
message("  UMAP done")


################################################################################
# STEP 9: Visualise
################################################################################

message("Generating UMAP plots ...")

make_plot <- function(grp, title, label = FALSE) {
  if (!grp %in% colnames(so@meta.data)) return(NULL)
  DimPlot(so, group.by = grp, reduction = "umap_harmony",
          label = label, repel = label, raster = TRUE) + ggtitle(title)
}

plots <- Filter(Negate(is.null), list(
  make_plot("seurat_clusters", "Clusters",      label = TRUE),
  make_plot("condition",       "By Condition"),
  make_plot("sample_id",       "By Sample"),
  make_plot("technology",      "By Technology"),
  make_plot("leiden_scVI",     "scVI Leiden")
))

n_cols <- min(2L, length(plots))
pdf(file.path(OUT_DIR, "BladderUrothelium_UMAP_overview.pdf"),
    width = 18, height = ceiling(length(plots) / n_cols) * 7)
print(wrap_plots(plots, ncol = n_cols))
dev.off()
message("  Saved: BladderUrothelium_UMAP_overview.pdf")

# Marker FeaturePlots — reload lognorm temporarily
markers <- c("Epcam", "Krt5", "Krt14", "Krt20", "Krt8", "Krt18",
             "Upk2", "Upk3a", "Trp63", "Upk1a", "Upk1b")
present <- intersect(markers, rownames(so))

if (length(present) > 0) {
  message("Restoring data layer for marker plots ...")
  so_tmp <- readRDS(IN_PATH)
  so[["RNA"]]$data <- so_tmp[["RNA"]]$data
  rm(so_tmp)
  gc()

  pdf(file.path(OUT_DIR, "BladderUrothelium_uro_markers.pdf"),
      width = 18, height = ceiling(length(present) / 5) * 6)
  print(FeaturePlot(so, features = present, reduction = "umap_harmony",
                    ncol = 5, raster = TRUE))
  dev.off()
  message("  Saved: BladderUrothelium_uro_markers.pdf")

  so[["RNA"]]$data <- NULL
  gc()
}


################################################################################
# STEP 10: Save
################################################################################

message(sprintf("Saving to %s ...", OUT_PATH))
saveRDS(so, OUT_PATH)

message("\n===== BladderUrothelium Harmony reclustering complete =====")
message(sprintf("Cells    : %s", format(ncol(so), big.mark = ",")))
message(sprintf("Clusters : %d", length(unique(so$seurat_clusters))))
message(sprintf("Output   : %s", OUT_PATH))
