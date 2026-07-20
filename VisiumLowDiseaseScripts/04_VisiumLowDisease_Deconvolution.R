################################################################################
# 04_VisiumLowDisease_Deconvolution.R
#
# RCTD deconvolution of the 10-sample disease/reference Visium object
# (VisiumLowDisease_harmony_integrated.rds: UUO, bilateral-IRI time course,
# 2 healthy adult 10x references -- see 01/02) using the adult
# RenalUrothelium reference (final_annotation_with_uro: mature nephron
# segments + Urothelium), same reference as
# SpatialScripts/02_VisiumHD_Deconvolution.R.
#
# Unlike VisiumLowDevelopmentScripts/05_VisiumLow_Deconvolution.R (which
# needed a stage-matched Chen2025 reference per sample because it spans
# embryo->aged), every sample here is adult mouse kidney -- a single global
# reference is appropriate, no per-stage subsetting needed.
#
# RenalUrothelium's RDS already carries a raw "counts" layer (unlike
# Chen2025's converted RDS, which needed a separate raw-count recovery step)
# -- confirmed integer-valued, 55,249 genes x 968,373 cells, RNA assay.
# 18,650 of the Visium object's 19,465 genes (the shared probe panel, see
# 01/02) overlap this reference; create.RCTD() intersects reference/query
# genes internally, no manual pre-filtering needed.
#
# spacexr::Reference() rejects "/" in cell type names (hit this exact issue
# with Chen2025's "B/Plasma" label in the dev-atlas script) -- defensively
# sanitized here too even though final_annotation_with_uro's labels
# (Asc-Vasa-Recta, B lymph, CD-Trans, ... Urothelium) don't currently
# contain one.
#
# Loops over all 10 image slots (each named "slice1.<sample_id>", except the
# first-merged sample which keeps the plain "slice1" -- see
# 02_VisiumLowDisease_Integrate_Harmony.R's merge() call), one RCTD run per
# sample, "full" mode (55um spots contain multiple cells).
#
# Input:  VisiumLowDiseaseScripts/output/VisiumLowDisease_harmony_integrated.rds
#         RenalUrotheliumScripts/output/RenalUrothelium_allcells_scvi_annotations_meta_uro_updated.rds
# Output: VisiumLowDisease_deconvolved.rds
#         VisiumLowDisease_DeconvCelltypeSpatialPlot.pdf
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(spacexr)
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(dplyr)
})

options(future.globals.maxSize = 16 * 1024^3)

log_mem <- function(label) {
  g  <- gc(verbose = FALSE)
  mb <- sum(g[, 2])
  message(sprintf("  [mem] %s: %.1f GB in use", label, mb / 1024))
}

SCRIPT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/VisiumLowDiseaseScripts"
OUT_DIR    <- file.path(SCRIPT_DIR, "output")
OBJECT_RDS <- file.path(OUT_DIR, "VisiumLowDisease_harmony_integrated.rds")
REF_RDS    <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/RenalUrotheliumScripts/output/RenalUrothelium_allcells_scvi_annotations_meta_uro_updated.rds"
OUT_DECONV <- file.path(OUT_DIR, "VisiumLowDisease_deconvolved.rds")
OUT_PDF    <- file.path(OUT_DIR, "VisiumLowDisease_DeconvCelltypeSpatialPlot.pdf")

# ── Load the disease/reference Visium object ────────────────────────────────
message("==> Loading VisiumLowDisease integrated object ...")
object <- readRDS(OBJECT_RDS)
log_mem("after loading Visium object")

# ── Build a single global RCTD reference from RenalUrothelium ──────────────
message("==> Building RCTD reference ...")
ref_so <- readRDS(REF_RDS)
log_mem("after loading ref_so")

ref_counts <- GetAssayData(ref_so, assay = "RNA", layer = "counts")
# spacexr::Reference() rejects "/" in cell type names -- see header comment.
cell_types_chr <- gsub("[/]", "_", as.character(ref_so$final_annotation_with_uro))
cell_types <- droplevels(as.factor(cell_types_chr))
names(cell_types) <- colnames(ref_so)
nUMI <- colSums(ref_counts)

rctd_ref <- spacexr::Reference(
  counts     = ref_counts,
  cell_types = cell_types,
  nUMI       = nUMI,
  min_UMI    = 1,
  n_max_cells = Inf
)
rm(ref_so, ref_counts, nUMI, cell_types); gc()
log_mem("after building RCTD reference")

message("Reference dimensions (genes x cells): ", paste(dim(rctd_ref@counts), collapse = " x "))
message("Cell type distribution:")
print(table(rctd_ref@cell_types))
ALL_CELLTYPES <- levels(rctd_ref@cell_types)

# ── Run RCTD per Visium sample against the shared reference ────────────────
message("\n==> Running RCTD per sample (full mode) ...")
weight_list <- vector("list", length(Images(object)))
names(weight_list) <- Images(object)

for (img_name in Images(object)) {
  slice_coords <- GetTissueCoordinates(object, image = img_name)
  cells <- rownames(slice_coords)
  samp  <- unique(object$sample_id[cells])
  stopifnot(length(samp) == 1)
  message(sprintf("  [%s] %s (%d spots)", img_name, samp, length(cells)))

  spatial_counts <- GetAssayData(object, assay = "Spatial", layer = "counts")[, cells]
  spatial_coords_df <- data.frame(
    x = slice_coords$imagecol,
    y = slice_coords$imagerow,
    row.names = cells
  )
  spatial_nUMI <- colSums(spatial_counts)

  query <- SpatialRNA(spatial_coords_df, spatial_counts, spatial_nUMI)
  myRCTD <- tryCatch(
    create.RCTD(query, rctd_ref, max_cores = 16, CELL_MIN_INSTANCE = 1),
    error = function(e) { message("    create.RCTD failed: ", e$message); NULL }
  )
  if (is.null(myRCTD)) next

  myRCTD <- run.RCTD(myRCTD, doublet_mode = "full")
  w <- as.matrix(normalize_weights(myRCTD@results$weights))

  # Reindex to the full reference cell-type schema (0-filled) so all 10
  # samples combine cleanly even if a given sample's fit happens to drop a
  # rare cell type column entirely.
  w_full <- matrix(0, nrow = nrow(w), ncol = length(ALL_CELLTYPES),
                    dimnames = list(rownames(w), ALL_CELLTYPES))
  w_full[, colnames(w)] <- w
  weight_list[[img_name]] <- w_full

  rm(myRCTD, query, spatial_counts); gc()
}
log_mem("after all RCTD runs")

# ── Combine per-sample weights, attach to the Seurat object ─────────────────
message("\n==> Combining RCTD weights across all samples ...")
rctd_weights <- do.call(rbind, weight_list[!sapply(weight_list, is.null)])
colnames(rctd_weights) <- paste0("rctd_", colnames(rctd_weights))

rctd_meta <- as.data.frame(rctd_weights)
rctd_meta$rctd_dominant_celltype <- sub(
  "^rctd_", "",
  colnames(rctd_weights)[max.col(rctd_weights, ties.method = "first")]
)

object <- AddMetaData(object, rctd_meta)
rm(weight_list, rctd_weights, rctd_meta); gc()
log_mem("after adding metadata")

# ── Save ─────────────────────────────────────────────────────────────────────
message("==> Saving deconvolved object ...")
saveRDS(object, OUT_DECONV, compress = FALSE)
message("  Saved: ", OUT_DECONV)

# ── Spatial plots of dominant cell type per sample ──────────────────────────
plot_spatial_manual <- function(obj, image_name, color_var, title = "",
                                 legend_title = color_var, pt_size = 1.2) {
  img_s4 <- obj@images[[image_name]]
  coords <- img_s4@coordinates
  coords$color_val <- obj[[color_var, drop = TRUE]][rownames(coords)]
  img_arr <- img_s4@image
  img_dim <- dim(img_arr)
  g <- rasterGrob(img_arr, width = unit(1, "npc"), height = unit(1, "npc"))

  ggplot(coords, aes(x = imagecol, y = imagerow)) +
    annotation_custom(g, xmin = 0, xmax = img_dim[2], ymin = 0, ymax = img_dim[1]) +
    geom_point(aes(color = color_val), size = pt_size) +
    scale_y_reverse() +
    coord_fixed(xlim = c(0, img_dim[2]), ylim = c(img_dim[1], 0)) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA),
          legend.text = element_text(size = 6)) +
    labs(color = legend_title) +
    ggtitle(title)
}

message("==> Plotting dominant cell type per sample ...")
pdf(OUT_PDF, width = 10, height = 8)
for (samp in unique(object$sample_id)) {
  message("  ", samp)
  sub_obj <- suppressWarnings(subset(object, cells = colnames(object)[object$sample_id == samp]))
  img_name <- Images(sub_obj)[1]
  tryCatch({
    print(
      plot_spatial_manual(sub_obj, img_name, "rctd_dominant_celltype",
                           title = samp, legend_title = "Dominant\ncell type")
    )
  }, error = function(e) {
    message("    Spatial plot failed for ", samp, ": ", e$message)
  })
}
dev.off()
message("  Saved: ", OUT_PDF)

message("\n==> Done.")
log_mem("final")
