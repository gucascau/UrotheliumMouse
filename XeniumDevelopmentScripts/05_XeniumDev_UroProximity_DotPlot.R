################################################################################
# 05_XeniumDev_UroProximity_DotPlot.R
#
# Two analyses on the 4-sample Xenium developmental-kidney dataset (GSE286051,
# W12/W92 x M/F), keyed off RCTD's per-cell dominant call
# (rctd_dominant_celltype) rather than a manually reconciled FinalCellType --
# unlike SpatialScripts' 12-sample UUO Xenium pipeline,
# 04_XeniumDev_MarkerHeatmap.R's cluster-annotation template has not been
# hand-filled yet, so RCTD's raw per-cell label is the only cell-type call
# available right now. Vocabulary differs from the UUO pipeline's
# FinalCellType: urothelium is "Urothelium" (not "Uro"), macrophages are
# "Macro" (not "Mac"), and there is no separate "Peri/VSMC" level -- only
# "Per" (pericyte; this reference doesn't split out VSMC). Confirmed against
# XeniumDev_RCTD_deconvolved.rds's dominant-celltype table
# (logs/xeniumdev_deconv_8293710.out).
#
# RCTD doublet-mode weights sum to 1 across at most 2 cell types per cell, and
# rctd_dominant_celltype is just argmax over them -- among cells called
# "Urothelium", median rctd_Urothelium weight is only 0.59 and 36.5% are
# below 0.5 (not even a majority, just a plurality). The dominant contaminant
# is Asc-Vasa-Recta, present in 37% of "Urothelium" cells (12,141 / 32,734) --
# biologically plausible, since urothelium lines the renal pelvis directly
# adjacent to vasa recta, so segmentation bleed between the two is expected
# in Xenium. Confirmed interactively (probe_rctd_purity.R). Every use of
# "Urothelium" below is therefore gated on rctd_Urothelium > MIN_URO_WEIGHT,
# not just rctd_dominant_celltype == "Urothelium" alone.
MIN_URO_WEIGHT <- 0.5

# 1) Proximity: nearest-neighbor centroid distance (Xenium tissue microns)
#    from each Urothelium cell to the closest cell of each nearby type, per
#    sample. Mirrors SpatialScripts/05_Xenium_Analyses.R's region-level
#    proximity block but at whole-sample scale. These 4 samples are
#    centroids-only FOVs (CreateCentroids + CreateFOV(type="centroids") in
#    02_XeniumDev_Integrate_Harmony.R -- the flat GEO deposit has no
#    cell_boundaries.parquet), so ImageDimPlot below renders Urothelium
#    centroids as points, not outlined/segmented cells like the UUO pipeline.
# 2) DotPlot: urothelium-cell marker expression across the 4 samples. This
#    Xenium panel is only 541 genes -- several requested markers (Krt14,
#    Krt15, Krt20, Cd24a, Spp1, Vcam1, Cd74, Fosl1) aren't in it.
#    Requested_MARKERS is intersected with rownames(object) and the drop
#    list is reported rather than erroring.
# 3) ImageDimPlot: per-sample spatial plot subset to Urothelium cells only.
# 4) DotPlot: 2 representative markers per rctd_dominant_celltype level,
#    across all cell types (not just Urothelium) -- see the celltype_markers
#    list below for the per-type gene picks and their caveats.
#
# Input:  XeniumDevelopmentScripts/output/XeniumDev_RCTD_deconvolved.rds
# Output: XeniumDev_UroProximity.csv
#         XeniumDev_UroProximitySummary.csv
#         XeniumDev_UroProximityBoxplot.pdf
#         XeniumDev_DotPlot_UroMarkers.pdf
#         XeniumDev_ImageDimPlot_UrotheliumOnly.pdf
#         XeniumDev_DotPlot_CellTypeMarkers.pdf
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

SCRIPT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/XeniumDevelopmentScripts"
OUT_DIR    <- file.path(SCRIPT_DIR, "output")
OBJECT_RDS <- file.path(OUT_DIR, "XeniumDev_RCTD_deconvolved.rds")

message("==> Loading RCTD-deconvolved Xenium object ...")
if (!file.exists(OBJECT_RDS)) {
  stop("Missing ", OBJECT_RDS, " -- run 03_XeniumDev_Deconvolution.R first.")
}
object <- readRDS(OBJECT_RDS)
DefaultAssay(object) <- "Xenium"

message("Dominant cell type distribution:")
print(table(object$rctd_dominant_celltype, useNA = "ifany"))

all_fovs <- Images(object)
message("Image slots found: ", paste(all_fovs, collapse = ", "))

# ── Proximity: Urothelium -> nearby cell types, per sample ─────────────────
proximity_targets <- c(Macrophage = "Macro", Fibroblast = "Fib",
                        Pericyte = "Per", TCell = "T lymph",
                        PrincipalCell = "PC", EndothelialCell = "Endo")

missing_targets <- setdiff(proximity_targets, unique(object$rctd_dominant_celltype))
if (length(missing_targets) > 0) {
  warning("Target cell type(s) not found in rctd_dominant_celltype: ",
          paste(missing_targets, collapse = ", "))
}

# Brute-force min distance (fine at these per-sample cell counts, low tens of
# thousands per type); returns NA for a source cell when the target type
# isn't present in that sample.
nearest_dist <- function(from_xy, to_xy) {
  if (nrow(from_xy) == 0) return(numeric(0))
  if (nrow(to_xy) == 0) return(rep(NA_real_, nrow(from_xy)))
  d <- sqrt(outer(from_xy[, 1], to_xy[, 1], "-")^2 + outer(from_xy[, 2], to_xy[, 2], "-")^2)
  apply(d, 1, min)
}

proximity_df <- do.call(rbind, lapply(all_fovs, function(fov_name) {
  coords <- GetTissueCoordinates(object, image = fov_name, which = "centroids")
  rownames(coords) <- coords$cell

  samp <- unique(object$sample_id[coords$cell])
  stopifnot(length(samp) == 1)
  message(sprintf("\n==> Sample: %s (fov %s, %d cells)", samp, fov_name, nrow(coords)))

  celltype   <- as.character(object$rctd_dominant_celltype[coords$cell])
  uro_weight <- object$rctd_Urothelium[coords$cell]
  uro_idx_all <- which(celltype == "Urothelium")
  uro_idx     <- which(celltype == "Urothelium" & uro_weight > MIN_URO_WEIGHT)
  if (length(uro_idx) == 0) {
    message("  [skip] no Urothelium cells above purity threshold")
    return(NULL)
  }
  uro_xy    <- as.matrix(coords[uro_idx, c("x", "y")])
  uro_cells <- coords$cell[uro_idx]
  message(sprintf("  %d Urothelium cells (%d before rctd_Urothelium > %.1f purity filter)",
                  length(uro_idx), length(uro_idx_all), MIN_URO_WEIGHT))

  do.call(rbind, lapply(names(proximity_targets), function(target_name) {
    target_type <- proximity_targets[[target_name]]
    target_idx  <- which(celltype == target_type)
    target_xy   <- as.matrix(coords[target_idx, c("x", "y")])
    data.frame(
      Sample     = samp,
      UroCell    = uro_cells,
      TargetType = target_name,
      nTarget    = length(target_idx),
      Distance   = nearest_dist(uro_xy, target_xy)
    )
  }))
}))

write.csv(proximity_df, file.path(OUT_DIR, "XeniumDev_UroProximity.csv"), row.names = FALSE)
message("\n  Saved proximity distance table.")

proximity_summary <- proximity_df %>%
  group_by(Sample, TargetType) %>%
  summarise(
    nUro       = sum(!is.na(Distance)),
    nTarget    = dplyr::first(nTarget),
    MedianDist = median(Distance, na.rm = TRUE),
    MeanDist   = mean(Distance, na.rm = TRUE),
    .groups    = "drop"
  )
write.csv(proximity_summary, file.path(OUT_DIR, "XeniumDev_UroProximitySummary.csv"), row.names = FALSE)
message("  Saved proximity summary table.")

sample_order <- sort(unique(proximity_df$Sample))
proximity_df$Sample <- factor(proximity_df$Sample, levels = sample_order)

ProximityBoxplot <- ggplot(proximity_df, aes(x = TargetType, y = Distance, fill = TargetType)) +
  geom_boxplot(outlier.size = 0.5, na.rm = TRUE) +
  facet_wrap(~ Sample, nrow = 2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(y = "Nearest-neighbor distance from Urothelium cell (µm)", x = NULL,
       title = "Distance from urothelial cells to nearby cell types")

ggsave(file.path(OUT_DIR, "XeniumDev_UroProximityBoxplot.pdf"), ProximityBoxplot, width = 7, height = 6)
message("  Saved proximity boxplot.")

# ── DotPlot: urothelium marker expression across samples ───────────────────
Requested_MARKERS <- c(
  # Uroplakins / urothelial keratins
  "Upk1a", "Upk1b", "Upk2", "Upk3a", "Upk3b",
  "Krt5", "Krt7", "Krt8", "Krt13", "Krt16", "Krt18", "Krt19", "Trp63",
  # Injury / inflammation
  "Havcr1", "Vcam1", "Cd74", "Fosl1",
  # Renal/urothelial developmental TFs and patterning
  "Pax2", "Pax8", "Six2", "Eya1", "Wt1", "Sox9",
  "Wnt3", "Sostdc1", "Shh", "Ptch1", "Bmp2", "Sfrp1", "Wif1",
  # Proliferation
  "Mki67", "Cenpf", "Stmn1", "Ube2c", "H2afx"
)
markers_present <- intersect(Requested_MARKERS, rownames(object))
markers_missing <- setdiff(Requested_MARKERS, rownames(object))
message(sprintf("\nMarkers present in 541-gene panel: %d / %d",
                length(markers_present), length(Requested_MARKERS)))
message("  Present: ", paste(markers_present, collapse = ", "))
if (length(markers_missing) > 0) {
  message("  Missing (not in panel): ", paste(markers_missing, collapse = ", "))
}

n_uro_all <- sum(object$rctd_dominant_celltype == "Urothelium", na.rm = TRUE)
UroOnly <- subset(object, subset = rctd_dominant_celltype == "Urothelium" & rctd_Urothelium > MIN_URO_WEIGHT)
UroOnly <- JoinLayers(UroOnly)
message(sprintf("Urothelium cells: %d (%d before rctd_Urothelium > %.1f purity filter)",
                ncol(UroOnly), n_uro_all, MIN_URO_WEIGHT))

if (length(markers_present) > 0 && ncol(UroOnly) > 0) {
  DotPlotUro <- DotPlot(
    UroOnly,
    features  = markers_present,
    group.by  = "sample_id",
    cols      = c("lightgrey", "red"),
    dot.scale = 8
  ) + RotatedAxis() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    guides(
      color = guide_colorbar(title = "AveExp"),
      size  = guide_legend(title  = "PerExp")
    )

  ggsave(file.path(OUT_DIR, "XeniumDev_DotPlot_UroMarkers.pdf"), DotPlotUro,
         width = max(6, 0.3 * length(markers_present) + 2), height = 4)
  message("  Saved DotPlot.")
}

# ── ImageDimPlot: Urothelium cells only, one page per sample ───────────────
# These FOVs are centroids-only (no segmentation polygons -- see header), so
# this is a point plot of Urothelium centroids over each sample's tissue
# extent, not an outlined-cell rendering.
message("\n==> Plotting Urothelium-only spatial distribution ...")
uro_fovs <- Images(UroOnly)
out_pdf_uro_spatial <- file.path(OUT_DIR, "XeniumDev_ImageDimPlot_UrotheliumOnly.pdf")
pdf(out_pdf_uro_spatial, width = 8, height = 7)
for (fov_name in uro_fovs) {
  fov_coords <- GetTissueCoordinates(UroOnly, image = fov_name, which = "centroids")
  samp <- unique(UroOnly$sample_id[fov_coords$cell])
  if (length(samp) == 0) next
  p <- tryCatch(
    ImageDimPlot(UroOnly,
      fov      = fov_name,
      group.by = "sample_id",
      cols     = "red",
      size     = 0.8,
      axes     = TRUE
    ) + ggtitle(paste("Urothelium cells only:", samp)) +
      theme(legend.position = "none"),
    error = function(e) { message("  Plot failed for ", samp, ": ", e$message); NULL }
  )
  if (!is.null(p)) print(p)
}
dev.off()
message("  Saved: ", out_pdf_uro_spatial)

# ── DotPlot: 2 representative markers per cell type (all identities) ───────
# Two genes per rctd_dominant_celltype level, picked from the 541-gene panel.
# Several nephron segments this panel can't tell apart get the same pair --
# LOH/CTAL/MTAL/MD all reduce to the TAL genes Umod+Kcnj1, since there's no
# cortex/medulla or macula-densa-specific probe in this panel -- and a few
# types (PEC, MC, DC, NK, the vasa recta subtypes) use the closest available
# proxy rather than a lineage-defining marker; both kinds of approximation
# are flagged inline below.
celltype_markers <- list(
  PTS1              = c("Lrp2", "Cubn"),
  PTS2              = c("Lrp2", "Slc22a6"),
  PTS3              = c("Lrp2", "Slc5a1"),
  PTS3T2            = c("Lrp2", "Slc22a8"),   # 2nd S3-family transporter -- panel can't separate further from PTS3
  DTL               = c("Aqp1", "Slc14a2"),
  ATL               = c("Cryab", "Cldn11"),
  "DTL-ATL"         = c("Aqp1", "Cryab"),     # transitional -- shares both thin-limb genes
  LOH               = c("Umod", "Kcnj1"),
  CTAL              = c("Umod", "Kcnj1"),     # same pair as LOH/MTAL -- no cortex/medulla discriminator in panel
  MTAL              = c("Umod", "Kcnj1"),
  MD                = c("Umod", "Kcnj1"),     # macula densa arises from TAL; no MD-specific probe in panel
  DCT               = c("Slc12a3", "Pvalb"),
  "DCT-CNT"         = c("Slc12a3", "Calb1"),
  CNT               = c("Calb1", "Scnn1b"),
  "CD-Trans"        = c("Aqp3", "Scnn1b"),
  PC                = c("Aqp2", "Avpr2"),
  ICA               = c("Slc4a1", "Atp6v1b1"),
  ICB               = c("Slc26a4", "Atp6v0d2"),
  Podo              = c("Nphs2", "Wt1"),
  PEC               = c("Cdh6", "Pax2"),      # approximate -- no PEC-specific probe in panel
  "Glom-Endo"       = c("Kdr", "Gpihbp1"),
  Endo              = c("Pecam1", "Tie1"),
  "Asc-Vasa-Recta"  = c("Plvap", "Gpihbp1"),
  "Desc-Vasa-Recta" = c("Aqp1", "Bst1"),
  "Vas-Afferens"    = c("Vwf", "Ackr3"),      # approximate -- no arterial/venous-specific probe in panel
  "Vas-Efferens"    = c("Vwf", "Eng"),
  Fib               = c("Col1a1", "Lum"),
  Per               = c("Rgs5", "Higd1b"),
  MC                = c("Des", "Itga8"),      # approximate -- overlaps Per
  Macro             = c("Cd68", "Marco"),
  DC                = c("Lamp3", "Clec4a3"),  # approximate
  Neutro            = c("Csf3r", "Cxcr2"),
  NK                = c("Klra8", "Ptprc"),    # approximate -- Ptprc (CD45) is pan-immune, not NK-specific
  "T lymph"         = c("Cd3e", "Cd3d"),
  "B lymph"         = c("Ms4a1", "Ptprc"),
  Urothelium        = c("Upk1b", "Krt5")
)

celltype_order           <- names(celltype_markers)
celltype_markers_vec     <- unique(unlist(celltype_markers, use.names = FALSE))
celltype_markers_present <- intersect(celltype_markers_vec, rownames(object))
message(sprintf("\nCell-type marker genes present: %d / %d",
                length(celltype_markers_present), length(celltype_markers_vec)))

object_typed <- subset(object, subset = !is.na(rctd_dominant_celltype))
Idents(object_typed) <- factor(object_typed$rctd_dominant_celltype, levels = rev(celltype_order))

DotPlotCellTypes <- DotPlot(
  object_typed,
  features  = celltype_markers_present,
  cols      = c("lightgrey", "blue"),
  dot.scale = 6
) + RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(
    color = guide_colorbar(title = "AveExp"),
    size  = guide_legend(title  = "PerExp")
  )

ggsave(file.path(OUT_DIR, "XeniumDev_DotPlot_CellTypeMarkers.pdf"), DotPlotCellTypes,
       width = max(10, 0.3 * length(celltype_markers_present) + 2), height = 10)
message("  Saved cell type marker DotPlot.")

message("\n==> Done.")
