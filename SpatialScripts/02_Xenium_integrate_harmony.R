################################################################################
# 02_Xenium_integrate_harmony.R
#
# Integrate 12 Xenium mouse kidney samples (UUO time-course) with Harmony:
#
#   GSM8325615  ShamL   — sham,   left  kidney
#   GSM8325616  ShamR   — sham,   right kidney
#   GSM8325617  Hour4L  — 4 h UUO,  left
#   GSM8325618  Hour4R  — 4 h UUO,  right
#   GSM8325619  Hour12L — 12 h UUO, left
#   GSM8325620  Hour12R — 12 h UUO, right
#   GSM8325621  Day2L   — day 2 UUO, left
#   GSM8325622  Day2R   — day 2 UUO, right
#   GSM8325623  Day14L  — day 14 UUO, left
#   GSM8325624  Day14R  — day 14 UUO, right
#   GSM8325625  Week6L  — week 6 UUO, left
#   GSM8325626  Week6R  — week 6 UUO, right
#
# Protocol reference:
#   10X Genomics Xenium 5K analysis journey (Seurat-based workflow)
#   https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2
#   https://colab.research.google.com/github/10XGenomics/analysis_guides/blob/main/Xenium_5k_data_analysis_journey.ipynb#scrollTo=V7gqfNcdKraN
#
# Steps:
#   1.  Extract *_output.tar.gz archives if not already done
#   2.  LoadXenium() per sample + QC filter (nCount_Xenium > 5, nFeature > 3)
#   3.  Merge all 12 samples
#   4.  SCTransform (assay = "Xenium")
#   5.  SketchData (LeverageScore, 5 000 cells per sample)
#   6.  RunPCA on sketch using all features
#   7.  RunHarmony on sketch (batch = sample_id)
#   8.  FindNeighbors + FindClusters + RunUMAP (harmony embedding)
#   9.  ProjectData back to full dataset
#  10.  Spatial + UMAP visualizations → PDF
#  11.  Save integrated object
################################################################################

library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)
library(dplyr)

options(future.globals.maxSize = 32 * 1024^3)

log_mem <- function(label) {
  g  <- gc(verbose = FALSE)
  mb <- sum(g[, 2])
  message(sprintf("  [mem] %s: %.1f GB in use", label, mb / 1024))
}

# ── Paths ─────────────────────────────────────────────────────────────────────
XEN_BASE <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSpatialData/Xenium"

OUT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/SpatialScripts/output"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OUT_RDS <- file.path(OUT_DIR, "Xenium_harmony_integrated.rds")
OUT_PDF <- file.path(OUT_DIR, "Xenium_harmony_clusters.pdf")

# ── Sample manifest ───────────────────────────────────────────────────────────
# inner_dir: the directory name packed inside each *_output.tar.gz
sample_meta <- data.frame(
  gsm       = c("GSM8325615","GSM8325616","GSM8325617","GSM8325618",
                "GSM8325619","GSM8325620","GSM8325621","GSM8325622",
                "GSM8325623","GSM8325624","GSM8325625","GSM8325626"),
  sample_id = c("ShamL","ShamR","Hour4L","Hour4R",
                "Hour12L","Hour12R","Day2L","Day2R",
                "Day14L","Day14R","Week6L","Week6R"),
  condition = c("Sham","Sham","Hour4","Hour4",
                "Hour12","Hour12","Day2","Day2",
                "Day14","Day14","Week6","Week6"),
  side      = c("L","R","L","R","L","R","L","R","L","R","L","R"),
  inner_dir = c(
    "ShamL-output-XETG00063__0005295__Region_5__20230717__191520",
    "ShamR-output-XETG00063__0005295__Region_6__20230717__191520",
    "Hour4L-output-XETG00063__0005161__Region_5__20230717__191520",
    "Hour4R-output-XETG00063__0005161__Region_6__20230717__191520",
    "Hour12L-output-XETG00063__0005295__Region_3__20230717__191520",
    "Hour12R-output-XETG00063__0005295__Region_4__20230717__191520",
    "Day2L-output-XETG00063__0005161__Region_3__20230717__191519",
    "Day2R-output-XETG00063__0005161__Region_4__20230717__191519",
    "Day14L-output-XETG00063__0005161__Region_1__20230717__191519",
    "Day14R-output-XETG00063__0005161__Region_2__20230717__191519",
    "Week6L-output-XETG00063__0005295__Region_1__20230717__191520",
    "Week6R-output-XETG00063__0005295__Region_2__20230717__191520"
  ),
  stringsAsFactors = FALSE
)

# ── Step 1: Extract archives ───────────────────────────────────────────────────
message("==> Extracting Xenium archives if needed ...")


# ── Step 2: Load each sample ──────────────────────────────────────────────────
message("==> Loading Xenium samples ...")

load_xenium_sample <- function(row) {
  data_dir <- file.path(XEN_BASE, row$inner_dir)
  message(sprintf("  Loading %s from %s ...", row$sample_id, data_dir))

  obj <- LoadXenium(data_dir, fov = "fov", segmentations = "cell")

  # QC filter: keep cells with detectable transcripts
  obj <- subset(obj, subset = nCount_Xenium > 5 & nFeature_Xenium > 3)

  # Add sample metadata
  obj$sample_id <- row$sample_id
  obj$condition <- row$condition
  obj$side      <- row$side
  obj$gsm       <- row$gsm

  log_mem(paste("after loading", row$sample_id))
  obj
}

sample_list <- vector("list", nrow(sample_meta))

sample_meta[1,]


for (i in seq_len(nrow(sample_meta))) {
  sample_list[[i]] <- load_xenium_sample(sample_meta[i, ])
}
names(sample_list) <- sample_meta$sample_id

# ── Step 3: Merge all samples ─────────────────────────────────────────────────
message("==> Merging all 12 samples ...")

object <- merge(
  x            = sample_list[[1]],
  y            = sample_list[-1],
  add.cell.ids = names(sample_list),
  project      = "Xenium_UUO_Integration"
)
rm(sample_list); gc()
log_mem("after merge")

# ── Step 4: SCTransform ───────────────────────────────────────────────────────
message("==> Running SCTransform ...")
# Xenium panels are pre-selected informative genes; normalize all features.
DefaultAssay(object) <- "Xenium"
object <- SCTransform(object, assay = "Xenium", clip.range = c(-10, 10),
                      verbose = TRUE)
log_mem("after SCTransform")

# ── Step 5: Geometric sketching ───────────────────────────────────────────────
message("==> Sketching data (LeverageScore, 5 000 cells per sample) ...")
object <- SketchData(
  object         = object,
  ncells         = 5000,
  method         = "LeverageScore",
  sketched.assay = "sketch"
)
DefaultAssay(object) <- "sketch"
log_mem("after sketching")

# ── Step 6: PCA on sketched data ──────────────────────────────────────────────
message("==> PCA on sketch ...")
# Use all features: Xenium panel is already curated, no HVG sub-selection needed
object <- ScaleData(object, assay = "sketch")
object <- RunPCA(object,
  assay          = "sketch",
  features       = rownames(object),
  reduction.name = "pca.sketch",
  npcs           = 30,
  verbose        = FALSE
)

# ── Step 7: Harmony on sketch (batch = sample_id) ────────────────────────────
message("==> Running Harmony on sketch ...")
object <- RunHarmony(
  object         = object,
  group.by.vars  = "sample_id",
  assay.use      = "sketch",
  reduction      = "pca.sketch",
  reduction.save = "harmony.sketch",
  theta          = 2,
  max_iter       = 20,
  verbose        = TRUE
)
log_mem("after Harmony")

# ── Step 8: Cluster + embed sketch ────────────────────────────────────────────
message("==> Clustering sketched cells ...")
object <- FindNeighbors(object,
  reduction = "harmony.sketch",
  dims      = 1:30,
  assay     = "sketch"
)
object <- FindClusters(object,
  resolution   = 0.3,
  cluster.name = "seurat_cluster.harmony.sketched"
)
object <- RunUMAP(object,
  reduction      = "harmony.sketch",
  dims           = 1:30,
  reduction.name = "umap.harmony.sketch",
  return.model   = TRUE
)
log_mem("after sketch clustering")

# ── Step 9: Project back to full dataset ──────────────────────────────────────
message("==> Projecting to full dataset ...")
object <- ProjectData(
  object             = object,
  assay              = "SCT",
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

# ── Step 10: Visualize ────────────────────────────────────────────────────────
message("==> Saving plots ...")
Idents(object) <- "seurat_cluster.harmony.projected"

pdf(OUT_PDF, width = 16, height = 7)

# UMAP: cluster identity
p_cluster <- DimPlot(object,
  reduction = "umap.harmony.sketch",
  group.by  = "seurat_cluster.harmony.projected",
  label     = TRUE, repel = TRUE, pt.size = 0.1
) + ggtitle("Xenium — Harmony clusters")

# UMAP: sample identity
p_sample <- DimPlot(object,
  reduction = "umap.harmony.sketch",
  group.by  = "sample_id",
  label     = FALSE, pt.size = 0.1
) + ggtitle("Sample identity")

print(p_cluster | p_sample)

# UMAP: condition
p_cond <- DimPlot(object,
  reduction = "umap.harmony.sketch",
  group.by  = "condition",
  label     = FALSE, pt.size = 0.1
) + ggtitle("Condition (UUO time-point)")

# UMAP: side (L = obstructed, R = contralateral)
p_side <- DimPlot(object,
  reduction = "umap.harmony.sketch",
  group.by  = "side",
  label     = FALSE, pt.size = 0.1
) + ggtitle("Side (L = obstructed, R = contralateral)")

print(p_cond | p_side)

# Spatial cluster plots for each sample
for (i in seq_len(nrow(sample_meta))) {
  samp <- sample_meta$sample_id[i]
  cells_samp <- WhichCells(object, expression = sample_id == samp)
  sub_obj <- subset(object, cells = cells_samp)
  Idents(sub_obj) <- "seurat_cluster.harmony.projected"
  tryCatch({
    print(
      ImageDimPlot(sub_obj, fov = "fov", cols = "polychrome",
                   axes = TRUE, size = 0.3) +
        ggtitle(paste("Spatial clusters:", samp))
    )
  }, error = function(e) {
    message("  ImageDimPlot failed for ", samp, ": ", e$message)
  })
}

dev.off()
message("  Saved: ", OUT_PDF)

# ── Step 11: Save ─────────────────────────────────────────────────────────────
message("==> Saving RDS ...")
saveRDS(object, OUT_RDS)
message("  Saved: ", OUT_RDS)

message("==> Done.")
log_mem("final")
