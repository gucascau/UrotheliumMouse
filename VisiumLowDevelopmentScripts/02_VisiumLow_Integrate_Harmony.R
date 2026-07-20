################################################################################
# 02_VisiumLow_Integrate_Harmony.R
#
# Integrate all 24 standard-resolution Visium developmental-kidney samples
# (GSE252772, see 01_VisiumLow_BuildMetadata.R) into one Harmony-corrected
# Seurat object, following this repo's established Visium integration
# pattern (SpatialScripts/01_VisiumHD_integrate_harmony_sketch.R): merge the
# raw "Spatial" counts assay (not the existing per-sample "SCT" assays, which
# were each fit independently and are not directly comparable across
# samples) -> NormalizeData -> FindVariableFeatures -> ScaleData -> RunPCA ->
# RunHarmony(batch = sample) -> cluster/UMAP.
#
# Unlike the VisiumHD integration, no SketchData/ProjectData step is used --
# total size here is 68,180 spots across 24 samples (vs. hundreds of
# thousands-to-millions of VisiumHD bins), small enough for Harmony directly
# on the full data.
#
# Metadata (Age, Sex, GSM, sample_id, mixed_sex_slide, n_sections_pooled) is
# attached to every spot from 01_VisiumLow_BuildMetadata.R's CSV before
# merging. Two caveats carried from that script:
#   - Sex = "Mixed" for the 5 P0 dual-sex slides (GSM8704846/849/850/851/852)
#     -- exclude these from any sex-stratified analysis.
#   - n_sections_pooled > 1 for the P0 dual-sex slides (2) and E16.5 slides
#     (4) -- those spots are not independent single-specimen replicates.
#
# Cell barcodes collide across samples (shared 10x barcode whitelist), so
# merge() uses add.cell.ids = sample_id.
#
# Input:  UsedSpatialData/VisiumLow/DevelopedKidney/*_obj.rds (24 files)
#         VisiumLowDevelopmentScripts/output/VisiumLow_sample_metadata.csv
#         (written by 01_VisiumLow_BuildMetadata.R)
# Output: VisiumLow_harmony_integrated.rds
#         VisiumLow_harmony_clusters.pdf
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(viridisLite)
})

options(future.globals.maxSize = 16 * 1024^3)

log_mem <- function(label) {
  g  <- gc(verbose = FALSE)
  mb <- sum(g[, 2])
  message(sprintf("  [mem] %s: %.1f GB in use", label, mb / 1024))
}

DATA_DIR   <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSpatialData/VisiumLow/DevelopedKidney"
SCRIPT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/VisiumLowDevelopmentScripts"
OUT_DIR    <- file.path(SCRIPT_DIR, "output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

STAGE_ORDER  <- c("E16.5", "P0", "W3", "W12", "W52", "W92")
STAGE_COLORS <- setNames(viridisLite::viridis(length(STAGE_ORDER)), STAGE_ORDER)

# ── Load metadata ────────────────────────────────────────────────────────────
meta_csv <- file.path(OUT_DIR, "VisiumLow_sample_metadata.csv")
if (!file.exists(meta_csv)) {
  stop("Missing ", meta_csv, " -- run 01_VisiumLow_BuildMetadata.R first.")
}
sample_meta <- read.csv(meta_csv, stringsAsFactors = FALSE)
sample_meta$Age <- factor(sample_meta$Age, levels = STAGE_ORDER)
message(sprintf("==> Loaded metadata for %d samples", nrow(sample_meta)))

# ── Load each sample, attach metadata, keep only the raw Spatial assay ──────
message("\n==> Loading 24 Visium objects ...")
obj_list <- vector("list", nrow(sample_meta))
names(obj_list) <- sample_meta$sample_id

for (i in seq_len(nrow(sample_meta))) {
  row <- sample_meta[i, ]
  message(sprintf("  [%d/%d] %s", i, nrow(sample_meta), row$sample_id))

  obj <- readRDS(file.path(DATA_DIR, row$rds_file))
  # These RDS files were saved with an older SeuratObject version whose
  # VisiumV1 image class lacks a "misc" slot the current version expects --
  # merge() fails validObject() on the raw image otherwise. Update first.
  obj <- suppressWarnings(UpdateSeuratObject(obj))
  # SCT is the default assay on load -- switch to Spatial before dropping SCT
  # (Seurat refuses to delete whichever assay is currently active). Not
  # comparable across independently SCTransformed samples anyway; we redo
  # normalization once on the merged object instead.
  DefaultAssay(obj) <- "Spatial"
  obj[["SCT"]] <- NULL

  # These objects store tissue_positions in full-resolution pixel space (e.g.
  # imagecol/imagerow up to ~18,000), but only the low-res image (~600x600 px)
  # is loaded into @image -- and unlike current Load10X_Spatial output, that
  # scaling isn't applied automatically at plot time for these older/updated
  # objects, so spots land far outside the visible image (confirmed by
  # comparing @image dims to @coordinates range, and by SpatialDimPlot
  # rendering spots invisible/off-frame). Pre-scale by scale.factors$lowres
  # once here so every downstream spatial plot is correctly aligned without
  # needing a per-plot patch.
  for (img_name in Images(obj)) {
    sf <- obj@images[[img_name]]@scale.factors$lowres
    obj@images[[img_name]]@coordinates$imagecol <- obj@images[[img_name]]@coordinates$imagecol * sf
    obj@images[[img_name]]@coordinates$imagerow <- obj@images[[img_name]]@coordinates$imagerow * sf
  }

  obj$GSM               <- row$GSM
  obj$sample_id         <- row$sample_id
  obj$Age               <- row$Age
  obj$Sex               <- row$Sex
  obj$mixed_sex_slide   <- row$mixed_sex_slide
  obj$n_sections_pooled <- row$n_sections_pooled

  obj_list[[i]] <- obj
}
log_mem("after loading all 24 samples")

# ── Merge (barcodes collide across samples -> add.cell.ids) ────────────────
message("\n==> Merging 24 samples ...")
object <- merge(
  x = obj_list[[1]],
  y = obj_list[-1],
  add.cell.ids = names(obj_list),
  project = "VisiumLow_DevelopedKidney"
)
rm(obj_list); gc()
log_mem("after merge")
message(sprintf("  Merged object: %d spots x %d genes across %d samples",
                ncol(object), nrow(object), length(unique(object$sample_id))))

# No JoinLayers step: these are old-style Assay (not Assay5) objects --
# merge() already concatenated the counts/data matrices directly, there are
# no per-sample "layers" to join (JoinLayers is Assay5-only and errors here).

# ── Normalize / HVG / scale on the merged raw counts ────────────────────────
message("\n==> Normalizing ...")
object <- NormalizeData(object, assay = "Spatial")

message("==> Finding variable features ...")
object <- FindVariableFeatures(object, assay = "Spatial", nfeatures = 3000)

message("==> Scaling data ...")
object <- ScaleData(object, assay = "Spatial")
log_mem("after ScaleData")

# ── PCA + Harmony (batch = sample) ──────────────────────────────────────────
message("\n==> PCA ...")
object <- RunPCA(object, assay = "Spatial", npcs = 30)

# Pinned to harmony 1.2.4 (pre-"Harmony2" rewrite; CRAN's harmony >= 2.0.0 is
# a rewrite with known stability issues and a different API) via
# remotes::install_version("harmony", version = "1.2.4") in the module R
# 4.4.0 personal library -- installed specifically for this reason, not the
# default harmony that was on the search path. Note even 1.2.4's
# RunHarmony.Seurat uses reduction.use/reduction.save, not the assay.use/
# reduction argument names in the VisiumHD precedent script
# (01_VisiumHD_integrate_harmony_sketch.R) -- that script's exact call was
# apparently never actually run against this environment's harmony version.
message("==> Running Harmony (batch = sample_id) ...")
object <- RunHarmony(
  object         = object,
  group.by.vars  = "sample_id",
  reduction.use  = "pca",
  reduction.save = "harmony",
  theta          = 2,
  max_iter       = 20,
  verbose        = TRUE
)
log_mem("after Harmony")

# ── Cluster + UMAP on the harmony embedding ─────────────────────────────────
message("\n==> Clustering ...")
object <- FindNeighbors(object, reduction = "harmony", dims = 1:30)
object <- FindClusters(object, resolution = 0.5, cluster.name = "seurat_cluster.harmony")
object <- RunUMAP(object, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
log_mem("after clustering")

message("\n  Spots per stage:")
print(table(object$Age))
message("  Spots per cluster:")
print(table(object$seurat_cluster.harmony))
message("  Cluster x stage crosstab:")
print(table(object$seurat_cluster.harmony, object$Age))

# ── Save ─────────────────────────────────────────────────────────────────────
# Diagnostic plots (UMAP + per-sample SpatialDimPlot x24) are a separate
# script (03_VisiumLow_Diagnostic_Plots.R) so they can be re-run/tweaked
# without repeating this ~10+ minute load/merge/Harmony/cluster step.
out_rds <- file.path(OUT_DIR, "VisiumLow_harmony_integrated.rds")
# compress = FALSE: gzip compression on an ~900 MB object (scale.data etc.)
# is slow; trade disk space for a much faster write.
saveRDS(object, out_rds, compress = FALSE)
message("\n  Saved: ", out_rds)

message("\n==> Done.")
