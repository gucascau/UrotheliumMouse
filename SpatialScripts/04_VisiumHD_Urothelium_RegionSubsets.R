################################################################################
# 04_VisiumHD_Urothelium_RegionSubsets.R
#
# Subset the RCTD-deconvolved VisiumHD object (both kidneys) into three
# manually-selected regions, save each as its own RDS, and plot urothelium
# biomarker spatial expression for each region.
#
# Region CSVs and their verified sample of origin (see
# resolve_3side_urothelium_region.py and inline notes below for how this was
# confirmed empirically — the barcode formats/sample assignment are NOT what
# an earlier half-finished exploration block in
# 01_VisiumHD_integrate_harmony_sketch.R assumed):
#
#   stromal–urothelial–myeloid .csv          -> kidney5p (slice1.008um.2)
#   Selected_UrotheliumRegionV2.csv           -> kidney5p (slice1.008um.2)
#   Selected_3SideUrotheliumRegion_V2.csv     -> kidney3p (slice1.008um)
#     (segmented-cell IDs, pre-resolved to square_008um bins via
#      Selected_3SideUrotheliumRegion_V2_resolved_008um_bins.txt)
#
# Cell name prefixes in the integrated object (from the merge() call in
# 01_VisiumHD_integrate_harmony_sketch.R, add.cell.ids = c("kidney3p","kidney")):
#   kidney3p sample -> "kidney3p_" prefix, image "slice1.008um"
#   kidney5p sample -> "kidney_"    prefix, image "slice1.008um.2"
#
# Input:  VisiumHD_kidney_both_deconvolved.rds
#         stromal–urothelial–myeloid .csv
#         Selected_UrotheliumRegionV2.csv
#         Selected_3SideUrotheliumRegion_V2_resolved_008um_bins.txt
# Output: VisiumHD_region_StromalUrothelialMyeloid.rds (+ PDF)
#         VisiumHD_region_UrotheliumRegionV2.rds (+ PDF)
#         VisiumHD_region_3SideUrothelium.rds (+ PDF)
################################################################################

library(Seurat)
library(ggplot2)
library(dplyr)

# ── Paths ─────────────────────────────────────────────────────────────────────
HD_BASE <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSpatialData/VisiumHD"
OUT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/SpatialScripts/output"

OBJECT_RDS <- file.path(OUT_DIR, "VisiumHD_harmony_KidneyOnlyintegrated.rds")

# Urothelial marker panel used throughout this project (see UrothelialMarkers
# in 01_VisiumHD_integrate_harmony_sketch.R line ~669)
UROTHELIAL_MARKERS <- c(
  "Krt5", "Krt14", "Krt20", "Trp63",
  "Upk1a", "Upk1b", "Upk2", "Upk3a", "Upk3b",
  "Foxa1", "Gata3", "Pparg",
  "Krt8", "Krt18", "Krt19"
)

# ── Load deconvolved object ────────────────────────────────────────────────────
message("==> Loading VisiumHD deconvolved object ...")
object <- readRDS(OBJECT_RDS)
# Object also carries a "sketch" assay (5000-bin/sample LeverageScore subsample
# used for fast Harmony clustering in 01_VisiumHD_integrate_harmony_sketch.R).
# DefaultAssay affects Cells(object)/subset() below — must be the full-resolution
# assay, or barcode matching silently checks against only the 10,000-bin sketch
# instead of the ~925,000 full-resolution bins (this was the actual cause of the
# "barcodes requested: 9957 | matched: 179"-type mismatches).
DefaultAssay(object) <- "Spatial.008um"
all_images <- Images(object)
message("Image slots found: ", paste(all_images, collapse = ", "))

markers_present <- intersect(UROTHELIAL_MARKERS, rownames(object))
message(sprintf("Urothelial markers present in object: %d / %d",
                length(markers_present), length(UROTHELIAL_MARKERS)))
message("  ", paste(markers_present, collapse = ", "))

# ── Helper: subset by barcode list, save RDS, plot markers ────────────────────
subset_region <- function(barcodes, prefix, image_name, region_name) {
  message(sprintf("\n==> Region: %s", region_name))

  cell_ids <- paste0(prefix, barcodes)
  matched  <- intersect(cell_ids, Cells(object))
  message(sprintf("  Barcodes requested: %d | matched in object: %d",
                  length(cell_ids), length(matched)))

  if (length(matched) == 0) {
    stop(sprintf(
      "No cells matched for region '%s' (prefix '%s', image '%s'). ",
      region_name, prefix, image_name),
      "Barcode prefix or sample assignment is likely wrong — do not proceed silently."
    )
  }

  region_obj <- subset(object, cells = matched)
  region_obj <- JoinLayers(region_obj)

  out_rds <- file.path(OUT_DIR, sprintf("VisiumHD_region_%s.rds", region_name))
  saveRDS(region_obj, out_rds)
  message("  Saved: ", out_rds)

  out_pdf <- file.path(OUT_DIR, sprintf("VisiumHD_region_%s_UrotheliumMarkers.pdf", region_name))
  pdf(out_pdf, width = 6, height = 5)
  for (gene in markers_present) {
    p <- SpatialFeaturePlot(region_obj,
      features       = gene,
      images         = image_name,
      crop           = TRUE,
      pt.size.factor = 12,
      image.alpha    = 1,
      alpha = c(0.1, 1)
    ) + ggtitle(paste(region_name, "-", gene)) +
      scale_fill_gradientn(
        colors   = c("transparent", "#FFFFD4", "#FED98E", "#FE8929", "#CC4C02"),
        na.value = "transparent"
      )
    print(p)
  }
  dev.off()

  message("  Saved: ", out_pdf)

  region_obj
}

# ── Region 1: the top urothelium region, that renal pelviscalyceal urotheliual line (kidney5p) ───────────────────────────
csv1 <- read.csv(file.path(HD_BASE, "Kidney1TopUroSection.csv"))
region1 <- subset_region(
  barcodes    = csv1$Barcode,
  prefix      = "kidney_",
  image_name  = "slice1.008um.2",
  region_name = "Pelviscalyceal_urothelium"
)
rm(region1); gc()

# ── Region 2: the middle urothelium region, that rrenal pelvis calyx lumen line (kidney5p) ────────────────────────
csv2 <- read.csv(file.path(HD_BASE, "Kidney1MiddleUroSection.csv"))
region2 <- subset_region(
  barcodes    = csv2$Barcode,
  prefix      = "kidney_",
  image_name  = "slice1.008um.2",
  region_name = "PelvisCalyxLumen_urothelium"
)
rm(region2); gc()


# ── Region 3: the bottom urothelium region, Suburothelial peri-urothelial region (kidney5p) ────────────────────────
csv3 <- read.csv(file.path(HD_BASE, "Kidney1MiddleUroSection.csv"))
region3 <- subset_region(
  barcodes    = csv3$Barcode,
  prefix      = "kidney_",
  image_name  = "slice1.008um.2",
  region_name = "Periurothelial_urothelium"
)
rm(region3); gc()


# ── Region 4: 3-side urothelium region (kidney3p, resolved from cell IDs) ─────
resolved_bins_file <- file.path(OUT_DIR, "Selected_3SideUrotheliumRegion_V2_resolved_008um_bins.txt")
bins4 <- readLines(resolved_bins_file)
region4 <- subset_region(
  barcodes    = bins4,
  prefix      = "kidney3p_",
  image_name  = "slice1.008um",
  region_name = "3SideUrothelium"
)
rm(region4); gc()

message("\n==> Done.")
