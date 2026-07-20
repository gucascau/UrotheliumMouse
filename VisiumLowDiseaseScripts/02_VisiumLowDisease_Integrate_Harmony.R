################################################################################
# 02_VisiumLowDisease_Integrate_Harmony.R
#
# Integrate the 10 adult-kidney disease-model / healthy-reference Visium
# samples (see 01_VisiumLowDisease_BuildMetadata.R) into one Harmony-corrected
# Seurat object. Deliberately kept separate from
# VisiumLowDevelopmentScripts' 24-sample developmental atlas (E16.5->W92, no
# injury) -- these 10 samples span UUO (GSE269063), a bilateral-IRI injury
# time course (GSE269884), and two healthy-adult 10x Genomics reference
# datasets, which is a different biological question from development.
#
# Three different on-disk formats, one per source (see 01's load_type):
#   - raw_geo_flat      : flat gzipped GEO supplementary files, no outs/
#                         folder -- gunzip into a staging dir shaped like a
#                         minimal spaceranger output, then Read10X() +
#                         Read10X_Image().
#   - spaceranger_outs  : a complete spaceranger outs/ dir already on disk --
#                         Read10X_h5() + Read10X_Image() directly, no staging.
#   - h5_plus_spatial_tar: filtered_feature_bc_matrix.h5 already on disk, but
#                         spatial/ is packed in a spatial.tar.gz -- untar just
#                         the spatial/ folder into a staging dir.
#
# image.type = "VisiumV1" is forced on every Read10X_Image() call. Seurat
# 5.5.1's default is "VisiumV2" (FOV/molecule-based, no @coordinates slot),
# but this repo's whole downstream spatial-plotting convention (manual
# @image + @coordinates ggplot rendering, established in
# VisiumLowDevelopmentScripts/03_VisiumLow_Diagnostic_Plots.R) depends on the
# legacy VisiumV1 slots -- confirmed both header (new-format
# tissue_positions.csv) and headerless (old-format tissue_positions_list.csv)
# inputs parse correctly under VisiumV1.
#
# Gene panel: Read10X_Image confirms (see conversation) the probe-based
# panel (19,465 genes; GSE269884 + Visium_FFPE_Mouse_Kidney) is an exact
# subset of the whole-transcriptome gene set (32,285 genes; GSE269063 +
# V1_Mouse_Kidney) -- same mm10-2020-A build throughout. Every sample is
# subset to the shared 19,465-gene probe panel right after loading, before
# merging, since Harmony cannot meaningfully batch-correct on genes that are
# structurally absent (not just undetected) from some samples.
#
# Unlike VisiumLowDevelopmentScripts (legacy Assay objects loaded from old
# RDS, no JoinLayers), every object here is freshly created with Seurat 5.5.1
# defaults -- Assay5 -- so merge() is followed by JoinLayers().
#
# Cell barcodes collide across samples (shared 10x barcode whitelist), so
# merge() uses add.cell.ids = sample_id.
#
# Input:  VisiumLowDiseaseScripts/output/VisiumLowDisease_sample_metadata.csv
#         (written by 01_VisiumLowDisease_BuildMetadata.R)
# Output: VisiumLowDisease_harmony_integrated.rds
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
})

options(future.globals.maxSize = 16 * 1024^3)

log_mem <- function(label) {
  g  <- gc(verbose = FALSE)
  mb <- sum(g[, 2])
  message(sprintf("  [mem] %s: %.1f GB in use", label, mb / 1024))
}

SCRIPT_DIR  <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/VisiumLowDiseaseScripts"
OUT_DIR     <- file.path(SCRIPT_DIR, "output")
STAGING_DIR <- file.path(OUT_DIR, "staging")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(STAGING_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Load metadata ────────────────────────────────────────────────────────────
meta_csv <- file.path(OUT_DIR, "VisiumLowDisease_sample_metadata.csv")
if (!file.exists(meta_csv)) {
  stop("Missing ", meta_csv, " -- run 01_VisiumLowDisease_BuildMetadata.R first.")
}
sample_meta <- read.csv(meta_csv, stringsAsFactors = FALSE)
message(sprintf("==> Loaded metadata for %d samples", nrow(sample_meta)))

# ── Per-load_type loaders -- all return a bare (unmerged) Seurat object ─────
load_raw_geo_flat <- function(row) {
  stage_dir <- file.path(STAGING_DIR, row$sample_id)
  spatial_dir <- file.path(stage_dir, "spatial")
  dir.create(spatial_dir, showWarnings = FALSE, recursive = TRUE)

  p <- row$matrix_prefix
  d <- row$matrix_dir
  file.copy(file.path(d, paste0(p, "_matrix.mtx.gz")),   file.path(stage_dir, "matrix.mtx.gz"), overwrite = TRUE)
  file.copy(file.path(d, paste0(p, "_barcodes.tsv.gz")), file.path(stage_dir, "barcodes.tsv.gz"), overwrite = TRUE)
  file.copy(file.path(d, paste0(p, "_features.tsv.gz")), file.path(stage_dir, "features.tsv.gz"), overwrite = TRUE)

  gunzip_to <- function(src, dest) {
    status <- system2("gunzip", c("-c", shQuote(src)), stdout = dest)
    if (status != 0) stop("gunzip failed for ", src)
  }
  gunzip_to(file.path(d, paste0(p, "_tissue_positions_list.csv.gz")), file.path(spatial_dir, "tissue_positions.csv"))
  gunzip_to(file.path(d, paste0(p, "_scalefactors_json.json.gz")),    file.path(spatial_dir, "scalefactors_json.json"))
  gunzip_to(file.path(d, paste0(p, "_tissue_lowres_image.png.gz")),   file.path(spatial_dir, "tissue_lowres_image.png"))
  gunzip_to(file.path(d, paste0(p, "_tissue_hires_image.png.gz")),    file.path(spatial_dir, "tissue_hires_image.png"))

  counts <- Read10X(data.dir = stage_dir)
  image  <- Read10X_Image(image.dir = spatial_dir, image.type = "VisiumV1")
  list(counts = counts, image = image)
}

load_spaceranger_outs <- function(row) {
  counts <- Read10X_h5(file.path(row$outs_dir, "filtered_feature_bc_matrix.h5"))
  image  <- Read10X_Image(image.dir = file.path(row$outs_dir, "spatial"), image.type = "VisiumV1")
  list(counts = counts, image = image)
}

load_h5_plus_spatial_tar <- function(row) {
  stage_dir <- file.path(STAGING_DIR, row$sample_id)
  dir.create(stage_dir, showWarnings = FALSE, recursive = TRUE)
  untar(row$spatial_tar, exdir = stage_dir)  # unpacks a top-level spatial/ dir

  counts <- Read10X_h5(row$h5_path)
  image  <- Read10X_Image(image.dir = file.path(stage_dir, "spatial"), image.type = "VisiumV1")
  list(counts = counts, image = image)
}

# ── Shared gene panel: the 19,465-gene probe set (see header comment) ──────
probe_panel_row <- sample_meta[sample_meta$n_genes_expected == 19465, ][1, ]
probe_panel_genes <- rownames(Read10X_h5(file.path(probe_panel_row$outs_dir, "filtered_feature_bc_matrix.h5")))
message(sprintf("==> Shared probe-panel gene set: %d genes", length(probe_panel_genes)))

# ── Load each sample, build a Seurat object, attach metadata, subset genes ──
message("\n==> Loading 10 Visium objects ...")
obj_list <- vector("list", nrow(sample_meta))
names(obj_list) <- sample_meta$sample_id

for (i in seq_len(nrow(sample_meta))) {
  row <- sample_meta[i, ]
  message(sprintf("  [%d/%d] %s (%s)", i, nrow(sample_meta), row$sample_id, row$load_type))

  loaded <- switch(row$load_type,
    raw_geo_flat        = load_raw_geo_flat(row),
    spaceranger_outs    = load_spaceranger_outs(row),
    h5_plus_spatial_tar = load_h5_plus_spatial_tar(row),
    stop("Unknown load_type: ", row$load_type)
  )
  counts <- loaded$counts
  image  <- loaded$image

  if (is.list(counts)) counts <- counts[["Gene Expression"]]
  stopifnot(nrow(counts) == row$n_genes_expected)

  obj <- CreateSeuratObject(counts = counts, assay = "Spatial", project = row$sample_id)

  # Full-resolution pixel coords loaded, but only the low-res image (~600px)
  # is attached -- same lowres-scale.factor pre-scaling fix documented in
  # VisiumLowDevelopmentScripts/02_VisiumLow_Integrate_Harmony.R, confirmed
  # necessary here too (imagecol/imagerow ranges far exceed the loaded
  # image's pixel dimensions before scaling).
  sf <- image@scale.factors$lowres
  image@coordinates$imagecol <- image@coordinates$imagecol * sf
  image@coordinates$imagerow <- image@coordinates$imagerow * sf

  image <- image[Cells(obj)]
  DefaultAssay(image) <- "Spatial"
  obj[["slice1"]] <- image

  # Subset to the shared probe-panel gene set before merging.
  shared_genes <- intersect(rownames(obj), probe_panel_genes)
  obj <- obj[shared_genes, ]

  obj$sample_id        <- row$sample_id
  obj$GSM               <- row$GSM
  obj$study             <- row$study
  obj$disease_model     <- row$disease_model
  obj$condition         <- row$condition
  obj$timepoint         <- row$timepoint
  obj$timepoint_hours   <- row$timepoint_hours
  obj$sex               <- row$sex
  obj$chemistry         <- row$chemistry
  obj$tissue_prep       <- row$tissue_prep

  obj_list[[i]] <- obj
}
log_mem("after loading all 10 samples")

# ── Merge (Assay5 -- JoinLayers needed; barcodes collide -> add.cell.ids) ──
message("\n==> Merging 10 samples ...")
object <- merge(
  x = obj_list[[1]],
  y = obj_list[-1],
  add.cell.ids = names(obj_list),
  project = "VisiumLowDisease"
)
rm(obj_list); gc()
object <- JoinLayers(object, assay = "Spatial")
log_mem("after merge + JoinLayers")
message(sprintf("  Merged object: %d spots x %d genes across %d samples",
                ncol(object), nrow(object), length(unique(object$sample_id))))

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

# Pinned to harmony 1.2.4, matching VisiumLowDevelopmentScripts/
# 02_VisiumLow_Integrate_Harmony.R -- see that script's comment for why
# (CRAN harmony >= 2.0.0 is an unstable "Harmony2" rewrite with a different
# API not installed here).
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

message("\n  Spots per study:")
print(table(object$study))
message("  Spots per sample:")
print(table(object$sample_id))
message("  Spots per cluster:")
print(table(object$seurat_cluster.harmony))
message("  Cluster x study crosstab:")
print(table(object$seurat_cluster.harmony, object$study))

# ── Save ─────────────────────────────────────────────────────────────────────
# Diagnostic plots are a separate script (03_VisiumLowDisease_Diagnostic_
# Plots.R) so they can be re-run/tweaked without repeating this load/merge/
# Harmony/cluster step.
out_rds <- file.path(OUT_DIR, "VisiumLowDisease_harmony_integrated.rds")
saveRDS(object, out_rds, compress = FALSE)
message("\n  Saved: ", out_rds)

unlink(STAGING_DIR, recursive = TRUE)

message("\n==> Done.")
