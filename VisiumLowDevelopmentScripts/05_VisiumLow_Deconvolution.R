################################################################################
# 05_VisiumLow_Deconvolution.R
#
# RCTD deconvolution of the 24-sample developmental Visium atlas
# (VisiumLow_harmony_integrated.rds), using the Chen2025 developmental
# atlas (MultiOmicSpatialMouseKidney_..._Chen2025NatGenet) as reference --
# NOT the adult RenalUrothelium reference used by
# SpatialScripts/02_VisiumHD_Deconvolution.R. That adult reference (mature
# nephron segments only: PTS1/2/3, CNT, DCT, etc.) has no embryonic/newborn
# progenitor states and would force E16.5/P0 spots into mismatched adult
# labels. Chen2025 is the exact single-cell/nucleus companion to this same
# Visium study (same animals, same 6 stages: E16.5, P0, W3, W12, W52, W92)
# and already carries per-stage-appropriate labels -- nephron progenitor
# (NP/NP_proliferate), ureteric bud progenitor (UBP), interstitial
# mesenchyme (IM/IM_proliferate) at E16.5/P0, mature segments expanding by
# W3+, and its own "Uro" (urothelial) label at every stage.
#
# Reference is built and matched PER STAGE, not once globally: for every
# Visium sample, RCTD runs against a Chen2025 reference subset to cells of
# that sample's own Age. This is the age-matching this repo's Chen2025-based
# scripts (UrotheliumDevelopmentScripts) already rely on -- see 01's header
# comment on why this Visium set and Chen2025 are treated as paired data.
#
# Chen2025's RDS only carries a log-normalized "data" layer (no counts
# survived its zellkonverter h5ad->Seurat conversion -- confirmed by
# inspecting @data directly: non-integer values, max ~8.3, consistent with
# log1p of a CPM-like normalization). RCTD needs raw UMI counts for its
# Poisson model. Raw counts were recovered from the original h5ad's
# CellxGene-schema `raw.X` slot (confirmed integer-valued, exact same
# cell/gene order as the converted RDS) via
# 00_export_Chen2025_rawcounts.py -- run that script first if
# output/Chen2025_rawcounts/ doesn't exist yet.
#
# Loops over all 24 image slots (slice1, slice1.2, ... slice1.24) in the
# integrated object, one RCTD run per sample against that sample's
# stage-matched reference, "full" mode (55um spots contain multiple cells).
# Per-sample weight matrices are reindexed to the full 26-level
# celltype_final schema (0-filled for cell types absent from that stage's
# reference -- see comment at combine step) before combining.
#
# Input:  VisiumLowDevelopmentScripts/output/VisiumLow_harmony_integrated.rds
#         VisiumLowDevelopmentScripts/output/Chen2025_rawcounts/ (from 00_*.py)
#         RawMouseSingleCellDatasets/..._Chen2025NatGenet_zellkonvertedConverted.rds
# Output: VisiumLow_deconvolved.rds
#         VisiumLow_DeconvCelltypeSpatialPlot.pdf
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

SCRIPT_DIR   <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/VisiumLowDevelopmentScripts"
OUT_DIR      <- file.path(SCRIPT_DIR, "output")
CHEN_COUNTS_DIR <- file.path(OUT_DIR, "Chen2025_rawcounts")
CHEN_RDS     <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/RawMouseSingleCellDatasets/MultiOmicSpatialMouseKidney_SexSpecific_lifespan_Chen2025NatGenet_zellkonvertedConverted.rds"
OBJECT_RDS   <- file.path(OUT_DIR, "VisiumLow_harmony_integrated.rds")
OUT_DECONV   <- file.path(OUT_DIR, "VisiumLow_deconvolved.rds")
OUT_PDF      <- file.path(OUT_DIR, "VisiumLow_DeconvCelltypeSpatialPlot.pdf")

STAGE_ORDER <- c("E16.5", "P0", "W3", "W12", "W52", "W92")

if (!dir.exists(CHEN_COUNTS_DIR)) {
  stop("Missing ", CHEN_COUNTS_DIR, " -- run 00_export_Chen2025_rawcounts.py first.")
}

# ── Load the Visium developmental atlas ─────────────────────────────────────
message("==> Loading VisiumLow integrated object ...")
object <- readRDS(OBJECT_RDS)
object$Age <- factor(object$Age, levels = STAGE_ORDER)
log_mem("after loading Visium object")

# ── Load Chen2025 reference: metadata from the RDS, raw counts from the ────
# ── exported h5ad raw.X (see header comment) ────────────────────────────────
message("==> Loading Chen2025 reference metadata + raw counts ...")
chen_meta_obj <- readRDS(CHEN_RDS)
chen_meta <- chen_meta_obj@meta.data[, c("Age", "celltype_final")]
rm(chen_meta_obj); gc()

chen_counts <- Read10X(data.dir = CHEN_COUNTS_DIR, unique.features = TRUE)
stopifnot(all(colnames(chen_counts) == rownames(chen_meta)))  # confirmed same order in 00_*.py
chen_meta$Age <- factor(chen_meta$Age, levels = STAGE_ORDER)
# spacexr::Reference() rejects "/" in cell type names (present in "B/Plasma")
# -- sanitize before building references, not just at the label that failed,
# since any future re-export could hit the same issue with other symbols.
chen_meta$celltype_final <- gsub("[/ ]", "_", as.character(chen_meta$celltype_final))
chen_meta$celltype_final <- droplevels(as.factor(chen_meta$celltype_final))
ALL_CELLTYPES <- levels(chen_meta$celltype_final)
log_mem("after loading Chen2025 reference")

message("  Chen2025 cells per stage:")
print(table(chen_meta$Age))

# ── Build one RCTD reference per developmental stage ────────────────────────
# Cell types with too few reference cells for a given stage must be dropped
# BEFORE spacexr::Reference() is built, not just via create.RCTD()'s
# CELL_MIN_INSTANCE argument -- discovered by debugging a "'dims' cannot be
# of length 0" failure (from RCTD's internal get_cell_mean()/sweep() choking
# on a 1-cell submatrix that loses its dim attribute) for every W12 and W52
# sample. Those two stages' Chen2025 references had LOH_AL_proliferating
# (n=1) and UBP (n=3) cells -- CELL_MIN_INSTANCE alone doesn't filter them
# out; process_cell_type_info requires the reference to already satisfy the
# minimum or it errors outright ("need a minimum of N cells for each cell
# type in the reference"). MIN_CELLS_PER_TYPE = 25 here, applied per stage.
#
# This does mean a rare-but-real identity (e.g. UBP, a nephric-progenitor-
# adjacent type) can be silently unavailable as an RCTD option in whichever
# stages it's too sparse to model reliably -- an intentional, unavoidable
# statistical limit (no reliable expression profile from a handful of
# cells), not a bug to paper over. Downstream analyses reading rctd_UBP (or
# any other column) should check which stages actually had enough reference
# cells for that type before comparing its "0" across stages, since a 0 here
# can mean either "genuinely absent" or "excluded, insufficient reference
# cells" -- these are not the same claim.
MIN_CELLS_PER_TYPE <- 25

message("\n==> Building per-stage RCTD references ...")
stage_refs <- list()
for (stage in STAGE_ORDER) {
  stage_cells_all <- rownames(chen_meta)[chen_meta$Age == stage]
  stage_types_all <- droplevels(chen_meta[stage_cells_all, "celltype_final"])
  type_counts <- table(stage_types_all)
  dropped_types <- names(type_counts)[type_counts < MIN_CELLS_PER_TYPE]
  if (length(dropped_types) > 0) {
    message(sprintf("  %s: dropping cell type(s) with < %d reference cells: %s",
                     stage, MIN_CELLS_PER_TYPE, paste(dropped_types, collapse = ", ")))
  }
  stage_cells <- stage_cells_all[!stage_types_all %in% dropped_types]
  message(sprintf("  %s: %d cells (%d after dropping rare types)",
                   stage, length(stage_cells_all), length(stage_cells)))

  stage_counts <- chen_counts[, stage_cells]
  stage_types  <- droplevels(chen_meta[stage_cells, "celltype_final"])
  names(stage_types) <- stage_cells
  stage_nUMI <- colSums(stage_counts)

  stage_refs[[stage]] <- spacexr::Reference(
    counts     = stage_counts,
    cell_types = stage_types,
    nUMI       = stage_nUMI,
    min_UMI    = 1,
    n_max_cells = Inf
  )
}
rm(chen_counts, chen_meta); gc()
log_mem("after building per-stage references")

# ── Run RCTD per Visium sample against its own stage-matched reference ─────
message("\n==> Running RCTD per sample (full mode) ...")
weight_list <- vector("list", length(Images(object)))
names(weight_list) <- Images(object)

for (img_name in Images(object)) {
  slice_coords <- GetTissueCoordinates(object, image = img_name)
  cells <- rownames(slice_coords)
  samp  <- unique(object$sample_id[cells])
  stage <- as.character(unique(object$Age[cells]))
  stopifnot(length(samp) == 1, length(stage) == 1)
  message(sprintf("  [%s] %s (stage %s, %d spots)", img_name, samp, stage, length(cells)))

  spatial_counts <- GetAssayData(object, assay = "Spatial", layer = "counts")[, cells]
  spatial_coords_df <- data.frame(
    x = slice_coords$imagecol,
    y = slice_coords$imagerow,
    row.names = cells
  )
  spatial_nUMI <- colSums(spatial_counts)

  query <- SpatialRNA(spatial_coords_df, spatial_counts, spatial_nUMI)
  myRCTD <- tryCatch(
    create.RCTD(query, stage_refs[[stage]], max_cores = 16, CELL_MIN_INSTANCE = MIN_CELLS_PER_TYPE),
    error = function(e) { message("    create.RCTD failed: ", e$message); NULL }
  )
  if (is.null(myRCTD)) next

  myRCTD <- run.RCTD(myRCTD, doublet_mode = "full")
  w <- as.matrix(normalize_weights(myRCTD@results$weights))

  # Reindex to the full 26-level celltype_final schema -- 0-filled for cell
  # types this stage's reference didn't include as an option (e.g. NP is
  # absent from the W92 reference, so W92 spots get 0, not NA -- they were
  # never offered that identity, which is the correct read for an argmax
  # dominant-cell-type call).
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

object <- readRDS(OUT_DECONV)
# ── Spatial plots of dominant cell type per sample ──────────────────────────
# SpatialDimPlot renders blank for these objects (see
# 03_VisiumLow_Diagnostic_Plots.R comment) -- manual ggplot rendering instead.
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
