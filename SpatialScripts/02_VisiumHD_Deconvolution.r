################################################################################
# 02_VisiumHD_Deconvolution.R
#
# RCTD deconvolution of VisiumHD 8 µm bins (both kidney samples) using
# annotated kidney scRNA-seq as reference (final_annotation_with_uro).
#
# Loops over every image slot in the KidneyOnly object (kidney3p + kidney5p)
# and runs RCTD per image, then combines weights across both samples.
#
# Input:  VisiumHD_harmony_KidneyOnlyintegrated.rds  (merged VisiumHD object)
#         RenalUrothelium_allcells_scvi_annotations_meta_uro_updated.rds
# Output: VisiumHD_kidney_deconvolved.rds (not this is for the 3p sample only)
#         DeconvCelltypeSpatialPlot.pdf
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

OBJECT_RDS  <- file.path(OUT_DIR, "VisiumHD_harmony_KidneyOnlyintegrated.rds")
REF_RDS     <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/RenalUrotheliumScripts/output/RenalUrothelium_allcells_scvi_annotations_meta_uro_updated.rds"
OUT_DECONV  <- file.path(OUT_DIR, "VisiumHD_kidney_deconvolved.rds")
OUT_PDF     <- file.path(OUT_DIR, "DeconvCelltypeSpatialPlot.pdf")

# ── Load integrated VisiumHD object ───────────────────────────────────────────
message("==> Loading VisiumHD integrated object ...")
object <- readRDS(OBJECT_RDS)
object <- JoinLayers(object, assay = "Spatial.008um")
log_mem("after loading object")

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

# ── Build SpatialRNA query from slice1.008um bins ─────────────────────────────
# Subset to kidney3p (slice1.008um) cells only — other samples have no image
# coordinates in this slot.
message("==> Building SpatialRNA query from slice1.008um ...")

slice_coords  <- GetTissueCoordinates(object, image = "slice1.008um", which = "centroids")
kidney3p_cells <- slice_coords$cell

spatial_counts <- GetAssayData(object, assay = "Spatial.008um", layer = "counts")
spatial_counts <- spatial_counts[, kidney3p_cells]

spatial_coords_df <- data.frame(
  x = slice_coords$x,
  y = slice_coords$y,
  row.names = slice_coords$cell
)
spatial_nUMI <- colSums(spatial_counts)

message(sprintf("  Query: %d genes x %d bins", nrow(spatial_counts), ncol(spatial_counts)))
query <- SpatialRNA(spatial_coords_df, spatial_counts, spatial_nUMI)
rm(slice_coords, spatial_counts, spatial_coords_df, spatial_nUMI); gc()
log_mem("after building query")

# ── Run RCTD ──────────────────────────────────────────────────────────────────
# "full" mode: 8 µm bins are smaller than a cell and can contain mixtures
message("==> Running RCTD (full mode) ...")
myRCTD <- create.RCTD(query, rctd_ref, max_cores = 16, CELL_MIN_INSTANCE = 1)
rm(query, rctd_ref); gc()

myRCTD <- run.RCTD(myRCTD, doublet_mode = "full")
log_mem("after RCTD run")

# ── Extract results and add to Seurat object ──────────────────────────────────
message("==> Adding RCTD weights to Seurat object ...")
rctd_weights <- as.matrix(normalize_weights(myRCTD@results$weights))
colnames(rctd_weights) <- paste0("rctd_", colnames(rctd_weights))

rctd_meta <- as.data.frame(rctd_weights)
rctd_meta$rctd_dominant_celltype <- sub(
  "^rctd_", "",
  colnames(rctd_weights)[max.col(rctd_weights, ties.method = "first")]
)

# AddMetaData matches on rownames — cells not in rctd_meta get NA automatically
object <- AddMetaData(object, rctd_meta)
rm(myRCTD, rctd_weights, rctd_meta); gc()
log_mem("after adding metadata")

# ── Save ──────────────────────────────────────────────────────────────────────
message("==> Saving deconvolved object ...")
saveRDS(object, OUT_DECONV)
message("  Saved: ", OUT_DECONV)

# ── Spatial plot of dominant cell type ────────────────────────────────────────
message("==> Plotting dominant cell type ...")
DeconvCelltypeSpatialPlot <- SpatialDimPlot(object,
  group.by       = "rctd_dominant_celltype",
  images         = "slice1.008um",
  crop           = TRUE,
  pt.size.factor = 2,
  image.alpha    = 1
) + theme(legend.position = "right")

ggsave(OUT_PDF, plot = DeconvCelltypeSpatialPlot, height = 6, width = 8)
message("  Saved: ", OUT_PDF)

message("==> Done.")
log_mem("final")
