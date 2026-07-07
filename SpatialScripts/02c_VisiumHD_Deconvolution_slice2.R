################################################################################
# 02c_VisiumHD_Deconvolution_slice2.R
#
# RCTD deconvolution of VisiumHD 8 µm bins for slice1.008um.2 (kidney5p)
# using annotated kidney scRNA-seq as reference (final_annotation_with_uro).
#
# Mirrors 02_VisiumHD_Deconvolution.r which targets slice1.008um (kidney3p).
#
# Input:  VisiumHD_harmony_KidneyOnlyintegrated.rds
#         RenalUrothelium_allcells_scvi_annotations_meta_uro_updated.rds
# Output: VisiumHD_kidney5p_deconvolved.rds
#         DeconvCelltypeSpatialPlot_slice2.pdf
################################################################################

library(Seurat)
library(spacexr)
library(ggplot2)
library(dplyr)

options(future.globals.maxSize = 16 * 1024^3)

log_mem <- function(label) {
  g  <- gc(verbose = FALSE)
  mb <- sum(g[, 2])
  message(sprintf("  [mem] %s: %.1f GB in use", label, mb / 1024))
}

# ── Paths ─────────────────────────────────────────────────────────────────────
OUT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/SpatialScripts/output"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OBJECT_RDS <- file.path(OUT_DIR, "VisiumHD_harmony_KidneyOnlyintegrated.rds")
REF_RDS    <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/RenalUrotheliumScripts/output/RenalUrothelium_allcells_scvi_annotations_meta_uro_updated.rds"
OUT_DECONV <- file.path(OUT_DIR, "VisiumHD_kidney5p_deconvolved.rds")
OUT_PDF    <- file.path(OUT_DIR, "DeconvCelltypeSpatialPlot_slice2.pdf")

TARGET_IMAGE <- "slice1.008um.2"

# ── Load VisiumHD integrated object ───────────────────────────────────────────
message("==> Loading VisiumHD integrated object ...")
object <- readRDS(OBJECT_RDS)
object <- JoinLayers(object, assay = "Spatial.008um")
log_mem("after loading object")

# Confirm target image exists
all_images <- Images(object)
message("Image slots found: ", paste(all_images, collapse = ", "))
if (!TARGET_IMAGE %in% all_images) {
  stop("Target image '", TARGET_IMAGE, "' not found. Available: ",
       paste(all_images, collapse = ", "))
}

# ── Build RCTD reference from annotated kidney scRNA-seq ──────────────────────
message("==> Building RCTD reference ...")
ref_so <- readRDS(REF_RDS)
log_mem("after loading ref_so")

ref_counts <- GetAssayData(ref_so, assay = "RNA", layer = "counts")
cell_types  <- droplevels(as.factor(ref_so$final_annotation_with_uro))
names(cell_types) <- colnames(ref_so)
nUMI <- colSums(ref_counts)

rctd_ref <- spacexr::Reference(
  counts     = ref_counts,
  cell_types = cell_types,
  nUMI       = nUMI,
  min_UMI    = 1,
  n_max_cells  = Inf
)
rm(ref_so, ref_counts, nUMI, cell_types); gc()
log_mem("after building RCTD reference")

message("Reference dimensions (genes x cells): ", paste(dim(rctd_ref@counts), collapse = " x "))
message("Cell type distribution:")
print(table(rctd_ref@cell_types))

# ── Build SpatialRNA query from slice1.008um.2 bins ───────────────────────────
message("==> Building SpatialRNA query from ", TARGET_IMAGE, " ...")

slice_coords  <- GetTissueCoordinates(object, image = TARGET_IMAGE, which = "centroids")
slice_cells   <- slice_coords$cell

spatial_counts <- GetAssayData(object, assay = "Spatial.008um", layer = "counts")
spatial_counts <- spatial_counts[, slice_cells]

spatial_coords_df <- data.frame(
  x = slice_coords$x,
  y = slice_coords$y,
  row.names = slice_cells
)
spatial_numi <- colSums(spatial_counts)

message(sprintf("  Query: %d genes x %d bins", nrow(spatial_counts), ncol(spatial_counts)))
query <- SpatialRNA(spatial_coords_df, spatial_counts, spatial_numi)
rm(slice_coords, spatial_counts, spatial_coords_df, spatial_numi); gc()
log_mem("after building query")

# ── Run RCTD ──────────────────────────────────────────────────────────────────
# "full" mode: 8 µm bins are smaller than a cell and can contain mixtures
message("==> Running RCTD (full mode) ...")
my_rctd <- create.RCTD(query, rctd_ref, max_cores = 20, CELL_MIN_INSTANCE = 1)
rm(query, rctd_ref); gc()

my_rctd <- run.RCTD(my_rctd, doublet_mode = "full")
log_mem("after RCTD run")

# ── Extract results and add to Seurat object ──────────────────────────────────
message("==> Adding RCTD weights to Seurat object ...")
rctd_weights <- as.matrix(normalize_weights(my_rctd@results$weights))
colnames(rctd_weights) <- paste0("rctd_", colnames(rctd_weights))

rctd_meta <- as.data.frame(rctd_weights)
rctd_meta$rctd_dominant_celltype <- sub(
  "^rctd_", "",
  colnames(rctd_weights)[max.col(rctd_weights, ties.method = "first")]
)

# AddMetaData matches on rownames; cells not in rctd_meta receive NA
object <- AddMetaData(object, rctd_meta)
rm(my_rctd, rctd_weights, rctd_meta); gc()
log_mem("after adding metadata")

message("Dominant cell type distribution:")
print(table(object$rctd_dominant_celltype, useNA = "ifany"))

# ── Save ──────────────────────────────────────────────────────────────────────
message("==> Saving deconvolved object ...")
saveRDS(object, OUT_DECONV)
message("  Saved: ", OUT_DECONV)

# ── Spatial plot ──────────────────────────────────────────────────────────────
message("==> Plotting dominant cell type ...")
deconv_plot <- SpatialDimPlot(object,
  group.by       = "rctd_dominant_celltype",
  images         = TARGET_IMAGE,
  crop           = TRUE,
  pt.size.factor = 2,
  image.alpha    = 1
) + ggtitle(paste("RCTD cell types:", TARGET_IMAGE)) +
  theme(legend.position = "right")

ggsave(OUT_PDF, plot = deconv_plot, height = 6, width = 8)
message("  Saved: ", OUT_PDF)

message("==> Done.")
log_mem("final")
