################################################################################
# 02b_VisiumHD_Deconvolution_BothKidneys.R
#
# RCTD deconvolution of VisiumHD 8 µm bins for BOTH kidney samples using
# annotated kidney scRNA-seq as reference (final_annotation_with_uro).
#
# Loops over every image slot in the KidneyOnly object (kidney3p + kidney5p)
# and runs RCTD per image, then combines weights across both samples.
#
# Strategy mirrors 03_Xenium_Deconvolution.R (per-sample loop).
#
# Input:  VisiumHD_harmony_KidneyOnlyintegrated.rds
#         RenalUrothelium_allcells_scvi_annotations_meta_uro_updated.rds
# Output: VisiumHD_kidney_both_deconvolved.rds
#         VisiumHD_RCTD_dominant_celltype.pdf
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
OUT_DECONV <- file.path(OUT_DIR, "VisiumHD_kidney_both_deconvolved.rds")
OUT_PDF    <- file.path(OUT_DIR, "VisiumHD_RCTD_dominant_celltype.pdf")

# ── Load VisiumHD integrated object ───────────────────────────────────────────
message("==> Loading VisiumHD integrated object ...")
object <- readRDS(OBJECT_RDS)
DefaultAssay(object) <- "Spatial.008um"
object <- JoinLayers(object, assay = "Spatial.008um")
log_mem("after loading object")

# Report all image slots — one per kidney sample in the KidneyOnly object
all_images <- Images(object)
message("Image slots found: ", paste(all_images, collapse = ", "))

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

# ── Run RCTD per image (one per kidney sample) ────────────────────────────────
# Load full count matrix once; subset per image inside the loop.
# "full" mode: 8 µm bins are smaller than a cell and can contain mixtures.
all_counts    <- GetAssayData(object, assay = "Spatial.008um", layer = "counts")
all_rctd_meta <- list()

for (img in all_images) {
  message(sprintf("\n==> RCTD for image: %s", img))

  img_coords <- GetTissueCoordinates(object, image = img, which = "centroids")
  img_cells  <- img_coords$cell

  coords_df <- data.frame(
    x = img_coords$x,
    y = img_coords$y,
    row.names = img_cells
  )

  img_counts <- all_counts[, img_cells]
  img_numi   <- colSums(img_counts)
  message(sprintf("  Bins: %d", length(img_cells)))

  query <- SpatialRNA(coords_df, img_counts, img_numi)
  rm(img_coords, coords_df, img_counts, img_numi); gc()

  my_rctd <- create.RCTD(query, rctd_ref, max_cores = 20, CELL_MIN_INSTANCE = 1)
  rm(query); gc()
  my_rctd <- run.RCTD(my_rctd, doublet_mode = "full")

  weights <- as.matrix(normalize_weights(my_rctd@results$weights))
  colnames(weights) <- paste0("rctd_", colnames(weights))

  meta_i <- as.data.frame(weights)
  meta_i$rctd_dominant_celltype <- sub(
    "^rctd_", "",
    colnames(weights)[max.col(weights, ties.method = "first")]
  )

  all_rctd_meta[[img]] <- meta_i
  rm(my_rctd, weights, meta_i); gc()
  log_mem(paste("after RCTD for", img))
}

rm(all_counts, rctd_ref); gc()

# Checkpoint before the combine step so RCTD results survive any downstream
# failure without needing to rerun the (many-hour) per-image RCTD loop.
saveRDS(all_rctd_meta, file.path(OUT_DIR, "all_rctd_meta_bothkidneys.rds"))

# ── Combine weights from both kidneys and add to Seurat object ────────────────
message("\n==> Adding RCTD weights to Seurat object ...")
# unname() prevents rbind.data.frame from prepending the list names (image
# IDs) onto each row name, which would double-prefix cell names (e.g.
# "slice1.008um.kidney3p_..." instead of "kidney3p_...") and make every row
# fail to match colnames(object) in AddMetaData below.
combined_meta <- do.call(rbind, unname(all_rctd_meta))
rm(all_rctd_meta); gc()

# AddMetaData matches on rownames; cells not in combined_meta receive NA
object <- AddMetaData(object, combined_meta)
rm(combined_meta); gc()
log_mem("after adding metadata")

message("Dominant cell type distribution:")
print(table(object$rctd_dominant_celltype, useNA = "ifany"))

# ── Save ──────────────────────────────────────────────────────────────────────
message("==> Saving deconvolved object ...")
saveRDS(object, OUT_DECONV)
message("  Saved: ", OUT_DECONV)

# ── Spatial plots — one panel per kidney sample ───────────────────────────────
message("==> Plotting ...")
pdf(OUT_PDF, width = 10, height = 7)

for (img in all_images) {
  tryCatch({
    p <- SpatialDimPlot(object,
      group.by       = "rctd_dominant_celltype",
      images         = img,
      crop           = TRUE,
      pt.size.factor = 2,
      image.alpha    = 1
    ) + ggtitle(paste("RCTD cell types:", img)) +
      theme(legend.position = "right")
    print(p)
  }, error = function(e) {
    message("  SpatialDimPlot failed for ", img, ": ", e$message)
  })
}

dev.off()
message("  Saved: ", OUT_PDF)

message("==> Done.")
log_mem("final")
