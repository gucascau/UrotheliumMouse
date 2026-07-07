################################################################################
# 01_VisiumHD_integrate_harmony_sketch.R
#
# Integrate three Visium HD mouse samples using Harmony + geometric sketching:
#   E-MTAB-15675-3            — K2a_322 bladder ctrl (16 µm bins)
#   Visium_HD_3prime_Mouse_Kidney — kidney 3′ capture (16 µm bins)
#   Visium_HD_Mouse_Kidney        — kidney 5′ capture (16 µm bins)
#
# Protocol: https://jef.works/blog/2025/04/22/harmony-with-sketching-in-seurat/
#
# Steps:
#   1. Extract tar.gz archives for the two kidney datasets (if needed)
#   2. Scaffold staging directory for E-MTAB-15675-3
#   3. Load all three with Load10X_Spatial(bin.size = 16)
#   4. Merge → NormalizeData → FindVariableFeatures → ScaleData
#   5. SketchData (LeverageScore, 5 000 cells per sample)
#   6. RunPCA on sketch
#   7. RunHarmony on sketch (batch = dataset)
#   8. FindNeighbors + FindClusters + RunUMAP (harmony embedding)
#   9. ProjectData back to full dataset
#  10. Save integrated object
################################################################################

library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)
library(dplyr)
library(jsonlite)
library(ggplot2)
# BiocManager::install("SpatialExperiment")
library(SpatialExperiment)
library(SummarizedExperiment)

library(spacexr)

options(future.globals.maxSize = 16 * 1024^3)

log_mem <- function(label) {
  g  <- gc(verbose = FALSE)
  mb <- sum(g[, 2])
  message(sprintf("  [mem] %s: %.1f GB in use", label, mb / 1024))
}

# ── Paths ─────────────────────────────────────────────────────────────────────
HD_BASE <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSpatialData/VisiumHD"

EMTAB_DIR   <- file.path(HD_BASE, "E-MTAB-15675-3")
KIDNEY3P_DIR <- file.path(HD_BASE, "Visium_HD_3prime_Mouse_Kidney")
KIDNEY_DIR   <- file.path(HD_BASE, "Visium_HD_Mouse_Kidney")

OUT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/SpatialScripts/output"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OUT_RDS  <- file.path(OUT_DIR, "VisiumHD_harmony_integrated.rds")
OUT_PDF  <- file.path(OUT_DIR, "VisiumHD_harmony_clusters.pdf")

# ── Step 1: Extract tar.gz archives ───────────────────────────────────────────
message("==> Extracting archives if needed ...")

# extract_if_needed <- function(dir, prefix) {
#   bin_dir <- file.path(dir, "binned_outputs", "square_016um")
#   if (!dir.exists(bin_dir)) {
#     message(sprintf("  Extracting %s ...", prefix))
#     tar_file <- file.path(dir, paste0(prefix, "_binned_outputs.tar.gz"))
#     if (!file.exists(tar_file)) stop("Cannot find: ", tar_file)
#     untar(tar_file, exdir = dir)

#     # Extract spatial tar (top-level images) if present
#     spatial_tar <- file.path(dir, paste0(prefix, "_spatial.tar.gz"))
#     if (file.exists(spatial_tar)) untar(spatial_tar, exdir = dir)

#     message("  Done.")
#   } else {
#     message(sprintf("  %s already extracted, skipping.", prefix))
#   }
# }

# extract_if_needed(KIDNEY3P_DIR, "Visium_HD_3prime_Mouse_Kidney")
# extract_if_needed(KIDNEY_DIR,   "Visium_HD_Mouse_Kidney")

# ── Step 2: Scaffold staging directory for E-MTAB-15675-3 ─────────────────────
# message("==> Setting up E-MTAB-15675-3 staging directory ...")

emtab_staging    <- file.path(EMTAB_DIR, "staging_008um")
bin008_dir       <- file.path(emtab_staging, "binned_outputs", "square_008um")
spatial_dir      <- file.path(bin008_dir, "spatial")
matrix_dir       <- file.path(bin008_dir, "filtered_feature_bc_matrix")
toplevel_spatial <- file.path(emtab_staging, "spatial")   # Load10X_Spatial reads images here

# setup_done_flag <- file.path(emtab_staging, ".setup_done")

# if (!file.exists(setup_done_flag)) {
#   dir.create(spatial_dir,      recursive = TRUE, showWarnings = FALSE)
#   dir.create(matrix_dir,       recursive = TRUE, showWarnings = FALSE)
#   dir.create(toplevel_spatial, recursive = TRUE, showWarnings = FALSE)

#   # Symlink matrix files
#   for (f in c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")) {
#     file.symlink(
#       from = normalizePath(file.path(EMTAB_DIR, f)),
#       to   = file.path(matrix_dir, f)
#     )
#   }

#   # Symlink images at top-level spatial/ (Load10X_Spatial reads images here)
#   for (img in c("tissue_hires_image.png", "tissue_lowres_image.png")) {
#     src <- file.path(EMTAB_DIR, paste0("K2a_322_16_bladder_ctrl_cropped_", img))
#     file.symlink(from = normalizePath(src), to = file.path(toplevel_spatial, img))
#   }

#   # Symlink positions + images inside binned spatial/ as well
#   file.symlink(
#     from = normalizePath(file.path(EMTAB_DIR,
#       "K2a_322_16_bladder_ctrl_cropped_tissue_positions.parquet")),
#     to   = file.path(spatial_dir, "tissue_positions.parquet")
#   )
#   for (img in c("tissue_hires_image.png", "tissue_lowres_image.png")) {
#     src <- file.path(EMTAB_DIR, paste0("K2a_322_16_bladder_ctrl_cropped_", img))
#     file.symlink(from = normalizePath(src), to = file.path(spatial_dir, img))
#   }

#   # scalefactors_json.json — standard Visium HD 08 µm defaults
#   scalefactors <- list(
#     spot_diameter_fullres     = 5.642,
#     tissue_hires_scalef       = 0.08181818181818182,
#     fiducial_diameter_fullres = 144,
#     tissue_lowres_scalef      = 0.024545454545454544
#   )
#   write_json(scalefactors,
#     file.path(spatial_dir, "scalefactors_json.json"),
#     auto_unbox = TRUE, pretty = TRUE)

#   file.create(setup_done_flag)
#   message("  Staging directory created.")
# } else {
#   message("  Staging directory already exists, skipping.")
# }

# ── Step 3: Load all three datasets ───────────────────────────────────────────
message("==> Loading Visium HD datasets at 16 µm bin size ...")

load_visiumhd <- function(data_dir, sample_name, bin = 8) {
  message(sprintf("  Loading %s ...", sample_name))


  obj <- Load10X_Spatial(data.dir = data_dir, bin.size = bin)
  obj$dataset     <- sample_name
  obj$sample_id   <- sample_name
  DefaultAssay(obj) <- paste0("Spatial.", sprintf("%03dum", bin))
  log_mem(paste("after loading", sample_name))
  obj
}
## Try Load10X_Spatial first; fall back to manual if it errors
obj_emtab <- tryCatch(
  load_visiumhd(emtab_staging, "E_MTAB_15675_3", bin = 8),
  error = function(e) {
    message("  Load10X_Spatial failed (", conditionMessage(e), "); using manual load ...")
    assay_name <- "Spatial.008um"
    slice_name <- "slice1.008um"
    counts <- Read10X(data.dir = matrix_dir)
    if (is.list(counts)) counts <- counts[["Gene Expression"]] %||% counts[[1]]
    image <- Read10X_Image(
      image.dir   = spatial_dir,
      image.name  = "tissue_hires_image.png",
      assay       = assay_name,
      slice       = slice_name,
      filter.matrix = TRUE
    )
    obj <- CreateSeuratObject(counts = counts, project = "Visium_HD_Bladder", assay = assay_name)
    image <- image[Cells(obj)]
    DefaultAssay(image) <- assay_name
    obj[[slice_name]] <- image
    DefaultAssay(obj) <- assay_name
    obj$dataset   <- "E_MTAB_15675_3"
    obj$sample_id <- "E_MTAB_15675_3"
    log_mem("after manual load E_MTAB_15675_3")
    obj
  }
)

# obj_emtab <- Load10X_Spatial(data.dir = emtab_staging, bin.size = 9) ## this is failed
obj_emtab$dataset   <- "Visium_HD_Bladder"
obj_emtab$sample_id <- "Visium_HD_Bladder"

#obj_emtab <- subset(obj_emtab, cells = which(obj_emtab$in_tissue == 1))

# Quick QC check
SpatialFeaturePlot(obj_emtab,
  features       = "nCount_Spatial.008um",
  crop           = TRUE,
  pt.size.factor = 1
)
library(jsonlite)

sf <- fromJSON(file.path(spatial_dir, "scalefactors_json.json"))
sf
sf$tissue_lowres_scalef
sf$tissue_hires_scalef

obj_emtab <- NormalizeData(obj_emtab)

SpatialFeaturePlot(obj_emtab, features = c("Upk2", "Upk1b"),
 crop = TRUE, pt.size.factor = 1) +
  theme(legend.position = "right")

# Check what image slots exist
names(obj_emtab@images)

# Check image dimensions (should be ~534 x 600 x 3 for a valid PNG)
img_slot <- names(obj_emtab@images)[1]
dim(obj_emtab@images[[img_slot]]@image)

# Check scale factors are set correctly
obj_emtab@images[[img_slot]]@scale.factors

# Just tissue + spots (no cluster coloring)
SpatialDimPlot(obj_emtab,
  crop           = TRUE,
  pt.size.factor = 1,
  label          = FALSE
)
# If you want to see only the tissue image with no spots at all:

SpatialDimPlot(obj_emtab,
  crop           = TRUE,
  pt.size.factor = 0,   # hides spots
  label          = FALSE
)
SpatialDimPlot(obj_emtab, crop = FALSE, pt.size.factor = 0)
obj_emtab@images[[img_slot]]@scale.factors
SpatialDimPlot(obj_emtab, crop = FALSE, pt.size.factor = 1.5)

SpatialFeaturePlot(obj_emtab,
  features       = c("Upk2", "Upk1b"),
  crop           = TRUE,
  pt.size.factor = 1,
  image.alpha    = 1
) + theme(legend.position = "right")


obj_kidney3p <- load_visiumhd(KIDNEY3P_DIR,   "Visium_HD_3prime_Kidney")
obj_kidney   <- load_visiumhd(KIDNEY_DIR,     "Visium_HD_Kidney")

# save them into RDS files for later use
saveRDS(obj_emtab, file = "obj_emtab.rds")
saveRDS(obj_kidney3p, file = "obj_kidney3p.rds")
saveRDS(obj_kidney, file = "obj_kidney.rds")

# ── Step 4: Merge → Normalize → Variable Features → Scale ─────────────────────
message("==> Merging datasets ...")
object <- merge(
  x = obj_emtab,
  y = list(obj_kidney3p, obj_kidney),
  add.cell.ids = c("bladder", "kidney3p", "kidney"),
  project = "VisiumHD_Integration"
)
rm(obj_emtab, obj_kidney3p, obj_kidney); gc()
log_mem("after merge")



ASSAY <- paste0("Spatial.", "008um")
DefaultAssay(object) <- ASSAY

message("==> Normalizing ...")
object <- NormalizeData(object, assay = ASSAY)

message("==> Finding variable features ...")
object <- FindVariableFeatures(object, nfeatures = 3000, assay = ASSAY)

message("==> Scaling data ...")
object <- ScaleData(object, assay = ASSAY)
log_mem("after ScaleData")

# ── Step 5: Geometric sketching ───────────────────────────────────────────────
message("==> Sketching data (LeverageScore, 5000 cells per dataset) ...")

# Sketch per dataset to preserve representation across batches
object <- SketchData(
  object        = object,
  ncells        = 5000,
  method        = "LeverageScore",
  sketched.assay = "sketch"
)

DefaultAssay(object) <- "sketch"
log_mem("after sketching")

# ── Step 6: PCA on sketched data ───────────────────────────────────────────────
message("==> PCA on sketch ...")
object <- FindVariableFeatures(object, assay = "sketch", nfeatures = 3000)
object <- ScaleData(object, assay = "sketch")
object <- RunPCA(object, assay = "sketch", reduction.name = "pca.sketch", npcs = 30)

# ── Step 7: Harmony batch correction on sketch ────────────────────────────────
message("==> Running Harmony on sketch (batch = dataset) ...")
object <- RunHarmony(
  object          = object,
  group.by.vars   = "dataset",
  assay.use       = "sketch",
  reduction       = "pca.sketch",
  reduction.save  = "harmony.sketch",
  theta           = 2,
  max_iter        = 20,
  verbose         = TRUE
)
log_mem("after Harmony")

# ── Step 8: Cluster and embed sketched cells (harmony embedding) ───────────────
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

# ── Step 9: Project back to full dataset ─────────────────────────────────────
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

# ── Step 10: Visualize and save ───────────────────────────────────────────────
message("==> Saving plots ...")
Idents(object) <- "seurat_cluster.harmony.projected"

pdf(OUT_PDF, width = 14, height = 6)

# UMAP colored by cluster
p_umap <- DimPlot(object,
  reduction = "umap.harmony.sketch",
  group.by  = "seurat_cluster.harmony.projected",
  label     = TRUE,
  repel     = TRUE
) + ggtitle("Harmony+Sketch clusters (UMAP)")

# UMAP colored by dataset
p_batch <- DimPlot(object,
  reduction = "umap.harmony.sketch",
  group.by  = "dataset",
  label     = FALSE
) + ggtitle("Dataset identity")

print(p_umap | p_batch)

# Spatial plots per sample
for (samp in c("Visium_HD_Bladder", "Visium_HD_3prime_Kidney", "Visium_HD_Kidney")) {
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


# check the meta data
object @meta.data %>% tail()

# the datasets from bladder is not good enough, now we will combine only the two kidneys. 

# ── Step 4: Merge → Normalize → Variable Features → Scale ─────────────────────
message("==> Merging datasets ...")

obj_kidney3p <- load_visiumhd(KIDNEY3P_DIR,   "Visium_HD_3prime_Kidney")
obj_kidney   <- load_visiumhd(KIDNEY_DIR,     "Visium_HD_Kidney")


object <- merge(
  x = obj_kidney3p,
  y = obj_kidney,
  add.cell.ids = c("kidney3p", "kidney"),
  project = "VisiumHD_KidneyIntegration"
)
rm(obj_kidney3p, obj_kidney); gc()
log_mem("after merge")



ASSAY <- paste0("Spatial.", "008um")
DefaultAssay(object) <- ASSAY

message("==> Normalizing ...")
object <- NormalizeData(object, assay = ASSAY)

message("==> Finding variable features ...")
object <- FindVariableFeatures(object, nfeatures = 3000, assay = ASSAY)

message("==> Scaling data ...")
object <- ScaleData(object, assay = ASSAY)
log_mem("after ScaleData")

# ── Step 5: Geometric sketching ───────────────────────────────────────────────
message("==> Sketching data (LeverageScore, 5000 cells per dataset) ...")

# Sketch per dataset to preserve representation across batches
object <- SketchData(
  object        = object,
  ncells        = 5000,
  method        = "LeverageScore",
  sketched.assay = "sketch"
)

DefaultAssay(object) <- "sketch"
log_mem("after sketching")

# ── Step 6: PCA on sketched data ───────────────────────────────────────────────
message("==> PCA on sketch ...")
object <- FindVariableFeatures(object, assay = "sketch", nfeatures = 3000)
object <- ScaleData(object, assay = "sketch")
object <- RunPCA(object, assay = "sketch", reduction.name = "pca.sketch", npcs = 30)


object@meta.data %>% head()
# ── Step 7: Harmony batch correction on sketch ────────────────────────────────
message("==> Running Harmony on sketch (batch = dataset) ...")
object <- RunHarmony(
  object          = object,
  group.by.vars   = "dataset",
  #assay.use       = "sketch",
  reduction       = "pca.sketch",
  reduction.save  = "harmony.sketch",
  theta           = 2,
  max_iter        = 20,
  verbose         = TRUE
)
log_mem("after Harmony")

# ── Step 8: Cluster and embed sketched cells (harmony embedding) ───────────────
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

# ── Step 9: Project back to full dataset ─────────────────────────────────────
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

# ── Step 10: Visualize and save ───────────────────────────────────────────────
message("==> Saving plots ...")
Idents(object) <- "seurat_cluster.harmony.projected"

pdf(OUT_PDF, width = 14, height = 6)

# UMAP colored by cluster
p_umap <- DimPlot(object,
  reduction = "umap.harmony.sketch",
  group.by  = "seurat_cluster.harmony.projected",
  label     = TRUE,
  repel     = TRUE
) + ggtitle("Harmony+Sketch clusters (UMAP)")

# UMAP colored by dataset
p_batch <- DimPlot(object,
  reduction = "umap.harmony.sketch",
  group.by  = "dataset",
  label     = FALSE
) + ggtitle("Dataset identity")

print(p_umap | p_batch)

# Spatial plots per sample
for (samp in c("Visium_HD_3prime_Kidney", "Visium_HD_Kidney")) {
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

message("==> Saving Kidney RDS ...")
saveRDS(object, "VisiumHD_harmony_KidneyOnlyintegrated.rds")
message("  Saved: ", OUT_RDS)

message("==> Kidney HD Done.")

message("==> Running Kidney cell deconvolution:")

object <- readRDS(file.path(OUT_DIR,"VisiumHD_harmony_KidneyOnlyintegrated.rds"))

all_images <- Images(object)

object <- JoinLayers(object, assay = "Spatial.008um")

# ── Build the RCTD reference from the annotated kidney scRNA allcells object ──
# celltype_final is the harmonized cell-type call, but only ~18% of cells have
# it filled in (the rest are NA) — RCTD needs a fully-labeled reference, so we
# restrict to the annotated subset.
KIDNEY_REF_RDS <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/RenalUrotheliumScripts/output/RenalUrothelium_allcells_scvi_annotations_meta_uro_updated.rds"

ref_so <- readRDS(KIDNEY_REF_RDS)
# annotation is final_annotation_with_uro
ref_counts   <- GetAssayData(ref_so, assay = "RNA", layer = "counts")
# ref_celltype <- setNames(droplevels(as.factor(ref_so$final_annotation_with_uro)), colnames(ref_so))
# Extract cell type annotations from the metadata (replace 'cell_type' with your actual metadata column name)
cell_types <- as.factor(ref_so$final_annotation_with_uro)
names(cell_types) <- colnames(ref_so)
# Ensure cell types are named vectors and drop any unused levels
cell_types <- droplevels(cell_types)

# Calculate nUMI (total counts per cell)
nUMI <- colSums(ref_counts)

# Build the RCTD Reference object
rctd_ref <- spacexr::Reference(counts = ref_counts, cell_types = cell_types, nUMI = nUMI,   min_UMI    = 1, n_max_cells  = Inf)

rm(ref_counts, nUMI, ref_so); gc()

# Examine reference. (optional)
print(dim(assay(rctd_ref))) # Gene expression matrix dimensions
#> [1] 750  75
table(colData(rctd_ref)$cell_types) # Number of occurrences of each cell type
#> 
dim(rctd_ref@counts)        # gene × cell dimensions
length(rctd_ref@cell_types) # number of cells
table(rctd_ref@cell_types)  # cells per type
head(rctd_ref@nUMI)         # nUMI per cellust@1120


# ── Build the SpatialRNA query from the non-urothelium kidney3p bins ─────────
# Restrict to the "slice1.008um" image (kidney3p sample), matching the scope
# of the urothelium correction above; the manually-called urothelium bins
# (common_cells) are excluded and keep their manual label instead of RCTD's.
slice_coords <- GetTissueCoordinates(object, image = "slice1.008um", which = "centroids")

spatial_counts <- GetAssayData(object, assay = "Spatial.008um", layer = "counts")

spatial_coords <- slice_coords[match(colnames(spatial_counts), slice_coords$cell), c("x", "y")]
rownames(spatial_coords) <- colnames(spatial_counts)
spatial_nUMI <- colSums(spatial_counts)

query <- SpatialRNA(spatial_coords, spatial_counts, spatial_nUMI)

# ── Run RCTD ──────────────────────────────────────────────────────────────────
# "full" mode: 8 um bins are far smaller than a cell, so each bin can be a
# mixture of many cell types rather than the 1-2 assumed by "doublet" mode.
myRCTD <- create.RCTD(query, rctd_ref, max_cores = 8, CELL_MIN_INSTANCE = 1)

myRCTD <- run.RCTD(myRCTD, doublet_mode = "full")

rctd_weights <- as.matrix(normalize_weights(myRCTD@results$weights))
colnames(rctd_weights) <- paste0("rctd_", colnames(rctd_weights))

obj_nonuro <- AddMetaData(obj_nonuro, as.data.frame(rctd_weights))
obj_nonuro$rctd_dominant_celltype <- sub(
  "^rctd_", "",
  colnames(rctd_weights)[max.col(rctd_weights, ties.method = "first")]
)

#saveRDS(myRCTD, file.path(OUT_DIR, "kidney3p_nonurothelium_RCTD.rds"))
saveRDS(object, file.path(OUT_DIR, "VisiumHD_kidney3p_urothelium_deconvolved.rds"))

DeconvCelltypeSpatialPlot <- SpatialDimPlot(object,
  group.by       = "deconv_celltype",
  images         = "slice1.008um",
  crop           = TRUE,
  pt.size.factor = 2,
  image.alpha    = 1
) + theme(legend.position = "right")

ggsave(file.path(OUT_DIR, "DeconvCelltypeSpatialPlot.pdf"),
       plot = DeconvCelltypeSpatialPlot, height = 6, width = 8)
message("==> Kidney cell deconvolution done. Saved: ", file.path(OUT_DIR, "VisiumHD_kidney3p_urothelium_deconvolved.rds"))


# We only concentrate on the UPK+ niches.
Idents(object) <- "seurat_cluster.harmony.projected"
# set up the order of the Identity
ObjectIdentityOrder<-  as.character(0:12)
# Reorder identities
levels(object) <- ObjectIdentityOrder


# we only show the identity 9
cells_ident9 <- WhichCells(object, idents = 9)

message("Images attached to object: ", paste(Images(object), collapse = ", "))

UrotheliumRegionHighlighSpatialHD <- SpatialDimPlot(object,
  cells.highlight = list(Urothelium = cells_ident9),
  cols.highlight   = c(Urothelium = "red", Unselected = "transparent"),
  pt.size.factor   = 3,
  image.alpha      = 1,
  alpha = 0.4,
  label            = F
) &
  scale_fill_manual(
    values = c(Urothelium = "red", Unselected = "transparent"),
    labels = c(Urothelium = "Urothelium", Unselected = "NonUrothelium"),
    name   = "Region"
  )

ggsave(file.path(OUT_DIR,"UrotheliumRegionHighlighSpatialHD.pdf"), plot = UrotheliumRegionHighlighSpatialHD, height = 8, width = 4)


# I have the selected regions that considered as the urothelium and draw the spatialFeature of these urothelium markers
# Urothelium markers
UrothelialMarkers<- c("Krt5", "Krt14", "Krt20",  "Trp63", "Upk1a","Upk1b", "Upk2", "Upk3a","Upk3b", "Foxa1", "Gata3", "Pparg", "Krt8", "Krt18", "Krt19")

p_UrothelialDotPlot_spatialHD <- DotPlot(object, features = UrothelialMarkers, cols = c("Spectral")) +
  ggtitle("Urothelium markers") +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(
    color = guide_colorbar(title = "AveExp"),
    size  = guide_legend(title = "PerExp")
  )
ggsave(paste0(OUT_DIR, "/UrothelialMarkersDotPlot_spatialHD.pdf"),
       plot = p_UrothelialDotPlot_spatialHD, width = 7, height = 3.5)

dir.create(file.path(OUT_DIR, "SpatialFeature"), showWarnings = FALSE, recursive = TRUE)
# Generete the spatial feature plot for the urothelium markers
for (gene in UrothelialMarkers) {
  #expr <- FetchData(object, vars = gene, layer = "data")[, 1]
  #object$plot_gene <- ifelse(expr > 0, expr, NA)

  p <- SpatialFeaturePlot(object, features = gene, pt.size.factor = 3, image.alpha = 1) +
  ggtitle(paste("Spatial Feature:", gene)) &
  scale_fill_gradientn(
    colors   = c("#FFFFD4", "#FED98E", "#FE8929", "#CC4C02"),
    na.value = "transparent"
  )
  ggsave(paste0(OUT_DIR, "/SpatialFeature/SpatialFeature_", gene, "_spatialHD.pdf"), plot = p, width = 6, height = 5)
}



### For the Renal cell type markers
mouse_kidney_markers <- list(
  Podocyte = c("Nphs1", "Nphs2", "Wt1", "Podxl"),
  Proximal_tubule = c("Slc34a1", "Lrp2", "Cubn", "Slc5a2", "Aqp1", "Hnf4a"),
  #Injured_repair_PT = c("Havcr1", "Lcn2", "Clu", "Spp1", "Sox9", "Vcam1", "Timp1", "Fn1"),
  Loop_of_Henle_TAL = c("Slc12a1", "Umod", "Kcnj1", "Bsnd"),
  DCT = c("Slc12a3", "Pvalb", "Trpm6", "Calb1"),
  CNT = c("Calb1", "Trpv5", "Slc8a1"),
  PrincipalCells = c("Aqp2", "Aqp3", "Aqp4", "Avpr2", "Fxyd4"),
  Intercalated = c("Foxi1", "Atp6v1b1", "Atp6v0d2", "Slc4a1", "Slc26a4"),
  Endothelial = c("Pecam1", "Cdh5", "Kdr", "Emcn", "Plvap"),
  Fibroblast = c("Col1a1", "Col1a2", "Col3a1", "Dcn", "Lum", "Pdgfra"),
  Pericyte_SMC = c("Rgs5", "Pdgfrb", "Acta2", "Tagln", "Myh11"),
  Immune = c("Ptprc", "Lyz2", "Adgre1", "Csf1r", "Cd3d", "Cd79a"),
)
mouse_kidney_markers_present <- lapply(
  mouse_kidney_markers,
  function(x) intersect(x, rownames(object))
)
mouse_kidney_markers_present<- c("Slc34a1", "Lrp2", "Cubn", "Slc12a1", "Umod", "Kcnj1","Slc12a3", "Pvalb", "Trpm6", "Calb1", "Trpv5", "Slc8a1","Aqp2", "Aqp3", "Aqp4","Atp6v1b1", "Atp6v0d2", "Slc4a1","Pecam1", "Cdh5", "Kdr", "Col1a2", "Col3a1", "Dcn","Acta2", "Tagln", "Myh11","Ptprc", "Lyz2", "Adgre1", "Csf1r", "Cd3d", "Cd79a")

p_RenalCellTypeMarkersDotPlot_spatialHD <- DotPlot(object, features = mouse_kidney_markers_present, cols = c("Spectral")) +
  ggtitle("Renal Cell Type markers") +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(
    color = guide_colorbar(title = "AveExp"),
    size  = guide_legend(title = "PerExp")
  )


ggsave(paste0(OUT_DIR, "/RenalCellTypeMarkersDotPlot_spatialHD.pdf"),
       plot = p_RenalCellTypeMarkersDotPlot_spatialHD, width = 9, height = 4)


### For the kidney injury markers
mouse_kidney_injury_genes <- c(
  "Spp1", "Havcr1", "Lcn2", "Clu", "Sox9", 
  "Vcam1", "Timp1", "Fn1", "Vim",
  "Krt6a", "Krt17",
  "Fos", "Jun", "Atf3"
)

p_InjuryMarkersDotPlot_spatialHD <- DotPlot(object, features = mouse_kidney_injury_genes, cols = c("Spectral")) +
  ggtitle("Injury markers") +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(
    color = guide_colorbar(title = "AveExp"),
    size  = guide_legend(title = "PerExp")
  )
dir.create(file.path(OUT_DIR, "SpatialFeature"), showWarnings = FALSE, recursive = TRUE)
# Generete the spatial feature plot for the urothelium markers

ggsave(paste0(OUT_DIR, "/RenalInjuryMarkersDotPlot_spatialHD.pdf"),
       plot = p_InjuryMarkersDotPlot_spatialHD, width = 6, height = 3.5)

for (gene in mouse_kidney_injury_genes) {
  p <- SpatialFeaturePlot(object, features = gene, pt.size.factor = 3, image.alpha = 1,   alpha = 0.4) +
    ggtitle(paste("Spatial Feature:", gene))
  ggsave(paste0(OUT_DIR, "/SpatialFeature/SpatialFeature_", gene, "_spatialHD.pdf"), plot = p, width = 6, height = 5)
}

# identify all the markers that within these regions
object <- JoinLayers(object)
UnBiasMarkersSpatialHD <- FindAllMarkers(object, test.use = "MAST", verbose = TRUE)

saveRDS(UnBiasMarkersSpatialHD, file ="UnBiasMarkersSpatialHD.rds")
UnBiasMarkersSpatialHD<- readRDS(file.path(OUT_DIR, "UnBiasMarkersSpatialHD.rds"))
# check the top 10 markers for each cluster, sorted by avg_log2FC, p adjust, and pct. expressed
top50_markers <- UnBiasMarkersSpatialHD %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), p_val_adj, desc(pct.1), .by_group = TRUE) %>%
  slice_head(n = 50) %>%
  ungroup() 

NovelMarkers<- UnBiasMarkersSpatialHD %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), p_val_adj, desc(pct.1), .by_group = TRUE) %>% filter(cluster == 9 & pct.1 >0.1 & pct.2 <0.05 & avg_log2FC > 0.5) %>% filter(!gene %in% UrothelialMarkers)%>% arrange(desc(pct.1), desc(avg_log2FC))  %>% slice_head(n = 15) %>% ungroup()


p_NovelBiomarkersDotPlot_spatialHD <- DotPlot(object, features = NovelMarkers$gene, cols = c("Spectral")) +
  ggtitle("Novel Candidate markers") +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(
    color = guide_colorbar(title = "AveExp"),
    size  = guide_legend(title = "PerExp")
  )
ggsave(paste0(OUT_DIR, "/RenalNovelMarkersDotPlot_spatialHD.pdf"),
       plot = p_NovelBiomarkersDotPlot_spatialHD, width = 6, height = 3.5)

# for the top novel markers
for (gene in NovelMarkers$gene) {
  expr <- GetAssayData(object, assay = "Spatial.008um", layer = "data")[gene, ]
object$plot_gene <- ifelse(expr > 0, expr, NA)

  p <- SpatialFeaturePlot(object, features = "plot_gene", pt.size.factor = 3, image.alpha = 1) +
  ggtitle(paste("Spatial Feature:", gene)) +
  scale_fill_gradientn(
    colors   = c("#FFFFD4", "#FED98E", "#FE8929", "#CC4C02"),
    na.value = "transparent"
  )
  ggsave(paste0(OUT_DIR, "/SpatialFeature/SpatialFeature_", gene, "_spatialHD.pdf"), plot = p, width = 6, height = 5)
}


# Upper-tract/renal-associated genes:
UpperRenalAssociatedGenes <- c("Pax8", "Pax2", "Glis3", "Fgfr2","Pkhd1", "Bicc1", "Magi1","Cgnl1", "Ptpn14","Col4a3", "Col4a4", "Col4a5")

p_UpperRenalAssociatedGenesDotPlot_spatialHD <- DotPlot(object, features = UpperRenalAssociatedGenes, cols = c("Spectral")) +
  ggtitle("Upper-Renal Associated Genes") +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(
    color = guide_colorbar(title = "AveExp"),
    size  = guide_legend(title = "PerExp")
  )
ggsave(paste0(OUT_DIR, "/RenalUpperRenalAssociatedGenesDotPlot_spatialHD.pdf"),
       plot = p_UpperRenalAssociatedGenesDotPlot_spatialHD, width = 6, height = 4)

# Functional enrichment analysis of DEGs from the cluster comparisions.



# Option 2: By x/y coordinate bounding box


coords <- GetTissueCoordinates(object, image = Images(object)[1], which = "centroids")

selected_cells <- coords$cell[
  coords$x > 3000 & coords$x < 8000 &
  coords$y > 5000 & coords$y < 12000
]

obj_region <- object
DefaultAssay(obj_region) <- "Spatial.008um"
obj_region@graphs <- list()
obj_region <- subset(obj_region, cells = selected_cells)
obj_region <- JoinLayers(obj_region)



# We used browser to identify the regions and further analyze thesem

# read the cells that contain these regions
UrotheliumRegionalCell1 <- read.csv(paste0(HD_BASE,"/", "Selected_UrotheliumRegionV2.csv"))
UrotheliumRegionalCell2 <- read.csv(paste0(HD_BASE,"/", "Selected_3SideUrotheliumRegion_V2.csv"))

# Get the coordnates of these cells
# section 1: 
# autromaticllly find the polygon image
Images(object)
poly_image2 <- grep("slice1.008um$", Images(object), value = TRUE)
# Get one x/y coordinate per segmented cell image 1
cell_coords2 <- GetTissueCoordinates(
  object = object,
  image = poly_image2,
  which = "centroids"
)
# Clean output
cell_coords2 <- as.data.frame(cell_coords2)

cell_coords2 %>% head()
# we need to add the string for "kidney3p_s_008um_"
UrotheliumRegionalCell2 %>% head()
UrotheliumRegionalCell1<-UrotheliumRegionalCell1 %>% mutate(Cell= paste0("kidney3p_",Barcode))
UrotheliumRegionalCell1 %>% nrow()
# UrotheliumRegionalCell2

#E standarize cell ID coluname 
# Clean output
cell_coords2 <- as.data.frame(cell_coords2)

cell_coords2 %>% head()
# Selected the common cells
common_cells <- intersect(UrotheliumRegionalCell1$Cell, cell_coords2$cell)
# 

head(cell_coords2$cell)           # actual IDs in the object
head(UrotheliumRegionalCell1$Cell) # IDs you constructed with "kidney3p_" prefix
length(common_cells)               # how many matched?


cell_coords2[common_cells,] %>% head()
common_cells_coordinates <- cell_coords2[cell_coords2$cell %in% common_cells, ]
common_cells_coordinates %>% head()
nrow(common_cells_coordinates)

# testcommoncells <- CreateSegmentation(common_cells_coordinates)

# object[["urothelium"]] <- Overlay(object[["slice1.008um"]], testcommoncells)
# uro_cells <- Cells(object[['urothelium']])

obj_tmp <- object
# DefaultAssay(obj_tmp) <- "Spatial.008um"   # switch away from sketch

# obj_tmp@graphs <- list()
urothelium <- subset(obj_tmp, cells = common_cells)
urothelium <- JoinLayers(urothelium)

SelectedMarkersSpatialFeaturePlots1 <- SpatialFeaturePlot(urothelium,
  features       = c("Upk2", "Upk1b", "Upk3b", "Upk3a"),
  crop           = TRUE,
  images         = "slice1.008um",
  pt.size.factor = 2,
  image.alpha    = 1
) + theme(legend.position = "right") + scale_fill_gradientn(
  colors   = c("transparent", "#FFFFD4", "#FED98E", "#FE8929", "#CC4C02"),
  na.value = "transparent"
)


#Fix 1: Make zero-expression spots transparent (best for seeing tissue)


SpatialFeaturePlot(urothelium,
  features       = "Upk3a",
  images         = "slice1.008um",
  crop           = TRUE,
  pt.size.factor = 2,
  image.alpha    = 1
) + scale_fill_gradientn(
  colors   = c("transparent", "#FFFFD4", "#FED98E", "#FE8929", "#CC4C02"),
  na.value = "transparent"
)
# Fix 2: Use min.cutoff to clip out low/zero values
SpatialFeaturePlot(urothelium,
  features       = "Upk3a",
  images         = "slice1.008um",
  crop           = TRUE,
  pt.size.factor = 1,
  image.alpha    = 1,
  min.cutoff     = "q10"   # bottom 10% (mostly zeros) gets lowest color
)

ggsave(file.path(OUT_DIR,"SelectedMarkersSpatialFeaturePlots_V1.pdf"), plot = SelectedMarkersSpatialFeaturePlots1, height = 15, width = 4)


# I want to perform a deconvolution for the spatial HD in the Kidney
# Here we manually corrected the urothleium cells but keep the other cell types

suppressPackageStartupMessages({
  library(spacexr)
})


# Option 1: From exported barcodes (Loupe Browser / browser tool)


# Load CSV of selected barcodes
selected <- read.csv("Selected_Region.csv")

# Build full cell IDs — match the format in your object
# Check format with: head(Cells(object))
selected_cells <- paste0("kidney3p_s_008um_", selected$Barcode)
selected_cells <- intersect(selected_cells, Cells(object))  # keep only valid ones

# Subset
obj_region <- object
DefaultAssay(obj_region) <- "Spatial.008um"
obj_region@graphs <- list()
obj_region <- subset(obj_region, cells = selected_cells)
obj_region <- JoinLayers(obj_region)

SpatialFeaturePlot(obj_region, features = c("Upk2", "Upk1b"), crop = TRUE, pt.size.factor = 4)
Option 2: By x/y coordinate bounding box


coords <- GetTissueCoordinates(object, image = Images(object)[1], which = "centroids")

selected_cells <- coords$cell[
  coords$x > 3000 & coords$x < 8000 &
  coords$y > 5000 & coords$y < 12000
]

obj_region <- object
DefaultAssay(obj_region) <- "Spatial.008um"
obj_region@graphs <- list()
obj_region <- subset(obj_region, cells = selected_cells)
obj_region <- JoinLayers(obj_region)
Option 3: Highlight selected cells without removing others (no subsetting)


# Keeps full tissue image but dims non-selected cells
SpatialDimPlot(object,
  cells.highlight = selected_cells,
  cols.highlight  = c("red", "grey90"),
  pt.size.factor  = 4,
  image.alpha     = 1
)
# For your case the pattern is always: build correct cell IDs → intersect to validate → clear graphs → subset → JoinLayers → plot. Option 3 is useful if you want to show context (where the region sits within the full tissue) without removing surrounding spots.

