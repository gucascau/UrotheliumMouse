################################################################################
# 05_Xenium_Analyses.R
#
# Subset the RCTD-deconvolved Xenium object into six manually-selected regions
# from Xenium Explorer, save each as its own RDS, and plot marker gene spatial
# expression for each region.  Follows the same structure as
# 04_VisiumHD_Urothelium_RegionSubsets.R.
#
# Region CSVs (Xenium Explorer exports, "Cell ID, Cluster, ..." columns):
#   Xenium_KidneyRight_Selection_{1,2,3}_coordinates.csv → ShamR, image fov.2
#   Xenium_KidneyCP_Selection_{1,2,3}_coordinates.csv    → currently empty;
#       fill in SAMPLE_CP / IMAGE_CP below once CSVs are re-exported.
#
# Cell name format in object: {sample_id}_{cell_id}
#   e.g.  ShamR_gplcidab-1
#
# Input:  output/Xenium_RCTD_deconvolved.rds
#         Xenium/SelectedRegions/*.csv
# Output: output/Xenium_region_{name}_{category}.rds  (+ PDF per region)
#         output/Xenium_UrotheliumSelectedRegions_{category}.rds
#         output/DotPlot_Xenium_UrotheliumSelectedRegions_{category}.pdf
#         output/DEGs_Xenium_UrotheliumSelectedRegions_{category}.csv
#         output/Heatmap_Top5DEGs_Xenium_UrotheliumSelectedRegions_{category}.pdf
################################################################################

library(Seurat)
library(ggplot2)
library(dplyr)
library(future)

# Crop()/Overlay() dispatch through future_lapply internally; the default
# future.globals.maxSize (500 MiB) is too small for this object's FOV/molecule
# data and silently falls back to the uncropped FOV. Raise it.
options(future.globals.maxSize = 4 * 1024^3)

# ── Paths ──────────────────────────────────────────────────────────────────────
XENIUM_BASE <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSpatialData/Xenium"
REGION_DIR  <- file.path(XENIUM_BASE, "SelectedRegions")
OUT_DIR     <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/SpatialScripts/output"

OBJECT_RDS <- file.path(OUT_DIR, "Xenium_RCTD_updated_deconvolved.rds")


# ── Marker gene panels ─────────────────────────────────────────────────────────
# Uncomment the block you want; the last active assignment wins.

RequestedCategoriesName <- "UrotheliumMarkers"
Requested_MARKERS <- c(
  "Upk1b", "Krt5", "Krt14","Krt15" ,"Krt19","Krt20"
)

RequestedCategoriesName <- "RenalInjuryMarkers"
Requested_MARKERS <- c(
  "Spp1", "Havcr1",  "Vcam1"
)

RequestedCategoriesName <- "NeighborCellsMarkers"
Requested_MARKERS <- c(
  "Aqp2",  # principal cells , "Scnn1g", 
  "Col1a2",  # Fibroblast "Col5a2",
  "Acta2", # VSMC , "Myh11",
  "Kdr", # Endo "Emcn"
  "Csf1r",# Mac ,"Il34",
  "Il2" # T cell "Il2ra" 
)

RequestedCategoriesName <- "EnrichedMarkers"
Requested_MARKERS <- c(
  "Gsdmc2", "Fxyd3", "Snx31", "Sprr1a", "Foxq1", "Psca", "S100a6", "Lgals3"
)

RequestedCategoriesName <- "UpperRenalUroGenes"
Requested_MARKERS <- c(
  "Pax8", "Pax2", "Glis3", "Fgfr1", "Pkhd1", "Bicc1",
  "Magi1", "Cgnl1", "Ptpn14", "Col4a3", "Col4a4", "Col4a5"
)

# ── Load deconvolved Xenium object ─────────────────────────────────────────────
message("==> Loading Xenium RCTD-deconvolved object ...")
object <- readRDS(OBJECT_RDS)

# check the annotation of the objects
object@meta.data %>% head()

DefaultAssay(object) <- "Xenium"
all_images <- Images(object)
message("Image slots found: ", paste(all_images, collapse = ", "))

markers_present <- intersect(Requested_MARKERS, rownames(object))
message(sprintf("Markers present in object: %d / %d",
                length(markers_present), length(Requested_MARKERS)))
message("  ", paste(markers_present, collapse = ", "))

# ── Fixed cell-type color map ───────────────────────────────────────────────
# Built once from every FinalCellType value in the full object (not per
# region/plot) so the same cell type always gets the same color across every
# ImageDimPlot page and every region — "cols = 'polychrome'" regenerates the
# palette per-call based only on the levels present in that call, which is
# why colors were drifting between plots.
#all_celltypes <- sort(unique(na.omit(object$FinalCellType)))

# Order by compartment (Epithelial -> Stromal -> Immune -> Unknown) instead of
# alphabetically, so the fixed legend/color order groups biologically related
# cell types together. Covers both the manual CellTypeAbbr vocabulary and the
# RCTD-only labels that can appear via the low-confidence rescue.
# epithelial_types <- c(
#   "PT-S1", "PT-S1/S2", "PT-S2/S3", "PT-S3", "PT-S3/TAL",
#   "FR-PT", "Inj-PT", "Inj-PT-S3", "Prolif-PT", "Prolif-Inj-PT",
#   "PTS2", "PTS3T2",
#   "TAL", "ATL", "DTL", "DTL-ATL",
#   "DCT", "DCT-CNT", "CNT",
#   "PC", "CD-Trans", "A-IC", "B-IC",
#   "MD", "Podo","Podo/GEnC", "PEC", "Uro"
# )
# stromal_types <- c(
#   "Fib", "Peri/VSMC", "Endo", "MC", "JG"
# )
# immune_types <- c("Mac", "Neutro", "T lymph", "B lymph", "NK", "DC")

# celltype_order <- c(epithelial_types, stromal_types, immune_types, "Unknown")

# # Anything in the data but not in the lists above falls at the end, flagged
# # so the lists can be updated rather than silently mis-ordered.
# uncategorized <- setdiff(all_celltypes, celltype_order)

all_celltypes <- c(  "PT-S1", "PT-S1/S2", "PT-S2/S3", "PT-S3", "PT-S3/TAL",
  "FR-PT", "Inj-PT", "Inj-PT-S3", "Prolif-PT", "Prolif-Inj-PT",
  "PTS2", "PTS3T2",
  "TAL", "ATL", "DTL", "DTL-ATL",
  "DCT", "DCT-CNT", "CNT",
  "PC", "CD-Trans", "A-IC", "B-IC",
  "MD", "Podo","Podo/GEnC", "PEC", "Uro",
  "Fib", "Peri/VSMC", "Endo", "MC", "JG",
  "Mac", "Neutro", "T lymph", "B lymph", "NK", "DC", "Unknown"
  )

celltype_colors <- setNames(
  scales::hue_pal()(length(all_celltypes)),
  all_celltypes
)


# ── Fixed marker (molecule) color map ───────────────────────────────────────
# Same problem, same fix, but for molecule dots: ImageDimPlot's molecule
# color palette (mols.cols) regenerates per call by default, so a given gene
# can render in a different color in the overview plot vs. its own per-gene
# page, and across different regions. Build one fixed gene->color map from
# Requested_MARKERS up front and reuse it everywhere.
marker_colors <- setNames(
  DiscretePalette(length(Requested_MARKERS), palette = "alphabet"),
  Requested_MARKERS
)

# Lock FinalCellType to a factor with the FULL level set (all cell types in
# the whole object), not just whichever subset happens to be non-empty in a
# given region. Some regions are missing certain cell types entirely; if
# FinalCellType stays a plain character vector, subset()/ImageDimPlot() infer
# levels per-plot from only what's present, and the name<->color pairing can
# drift for absent-elsewhere types. A fixed-level factor is carried through
# subset(), JoinLayers(), and Crop() unchanged, so every region plot uses the
# exact same name-to-color assignment even when a level has zero cells there.
object$FinalCellType <- factor(object$FinalCellType, levels = all_celltypes)

# ── Sample ↔ image mapping for selected regions ────────────────────────────────
# KidneyRight  → ShamR sample, confirmed by barcode matching (1638 / 1642 hit)
# KidneyLeft     → files currently empty; update prefix/fov when CSVs are ready
SAMPLE_RIGHT <- "ShamR_"
IMAGE_RIGHT  <- "fov.2"

SAMPLE_LEFT    <- "ShamL_"   # TODO: update to the correct sample_id
IMAGE_LEFT     <- "fov"    # TODO: update to the correct fov

object@meta.data %>% head()
# ── Helper: read Xenium Explorer CSV (skips # comment header lines) ───────────
read_xenium_csv <- function(csv_path) {
  info <- file.info(csv_path)
  if (is.na(info$size) || info$size == 0) {
    message("  [skip] Empty or missing: ", basename(csv_path))
    return(NULL)
  }
  df <- read.csv(csv_path, comment.char = "#", header = TRUE,
                 check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(df) == 0 || !"Cell ID" %in% colnames(df)) {
    message("  [skip] No 'Cell ID' column in: ", basename(csv_path))
    return(NULL)
  }
  df[["Cell ID"]]
}

# ── Helper: subset region, crop FOV, save RDS, plot segmentation + molecules ───
subset_region <- function(cell_ids, sample_prefix, image_name, region_name) {

  # cell_ids <- ids_kr1
  # sample_prefix <- SAMPLE_RIGHT
  # image_name <- IMAGE_RIGHT
  # region_name <- "KidneyRight_Sel1"

  message(sprintf("\n==> Region: %s", region_name))
  if (is.null(cell_ids)) {
    message("  [skip] No cell IDs — CSV was empty.")
    return(invisible(NULL))
  }

  full_ids <- paste0(sample_prefix, cell_ids)
  matched  <- intersect(full_ids, Cells(object))
  message(sprintf("  Cell IDs requested: %d | matched in object: %d",
                  length(full_ids), length(matched)))

  if (length(matched) == 0) {
    warning(sprintf(
      "No cells matched for region '%s' (prefix '%s', image '%s'). ",
      region_name, sample_prefix, image_name),
      "Check that the sample_prefix and image_name are correct.")
    return(invisible(NULL))
  }

  region_obj <- subset(object, cells = matched)
  region_obj <- JoinLayers(region_obj)

  out_rds <- file.path(OUT_DIR, sprintf("Xenium_region_%s_%s.rds",
                                         region_name, RequestedCategoriesName))
  saveRDS(region_obj, out_rds)
  message("  Saved: ", out_rds)

  # Crop the FOV to the bounding box of the selected cells (+ 50 µm padding),
  # then set DefaultBoundary to "segmentation" to show cell outlines.
  coords <- GetTissueCoordinates(region_obj, image = image_name, which = "centroids")
  pad   <- 10
  x_rng <- range(coords$x) + c(-pad, pad)
  y_rng <- range(coords$y) + c(-pad, pad)
  message(sprintf("  Crop bounds: x=[%.0f, %.0f]  y=[%.0f, %.0f]",
                  x_rng[1], x_rng[2], y_rng[1], y_rng[2]))

  cropped <- tryCatch(
    Crop(region_obj[[image_name]], x = x_rng, y = y_rng, coords = "tissue"),
    error = function(e) {
      message("  Crop() failed: ", e$message, " — using uncropped FOV")
      region_obj[[image_name]]
    }
  )
  region_obj[["zoom"]] <- cropped
  DefaultBoundary(region_obj[["zoom"]]) <- "segmentation"
  # check the molecular exists
  # # original image (before crop) — should show molecules present
  # Molecules(region_obj[[image_name]])
  # length(region_obj[[image_name]][["molecules"]])

  # # cropped "zoom" FOV — what ImageDimPlot(fov = "zoom", ...) actually uses
  # Molecules(region_obj[["zoom"]])
  # length(region_obj[["zoom"]][["molecules"]])
  # zoom_names <- names(region_obj[["zoom"]][["molecules"]])
  # setdiff(markers_present, zoom_names)


  out_pdf <- file.path(OUT_DIR, sprintf("Xenium_region_%s_%s.pdf",
                                         region_name, RequestedCategoriesName))
  pdf(out_pdf, width = 8, height = 6)

  # Page 1: overview — cell-type segmentation + all marker molecule dots
  p_all <- tryCatch(
    ImageDimPlot(region_obj,
      fov          = "zoom",
      group.by     = "FinalCellType",
      molecules    = markers_present,
      nmols        = 20000,
      mols.size    = 0.3,
      mols.cols    = marker_colors[markers_present],
      border.color = "white",
      border.size  = 0.1,
      cols         = celltype_colors,
      coord.fixed  = FALSE,
      axes         = TRUE
    ) + ggtitle(paste(region_name, "— cell types + all markers")) +
      theme(legend.position = "right"),
    error = function(e) { message("  Overview plot failed: ", e$message); NULL }
  )
  if (!is.null(p_all)) print(p_all)

  # One page per gene: cell-type segmentation + that gene's molecule dots
  for (gene in markers_present) {
    p <- tryCatch(
      ImageDimPlot(region_obj,
        fov          = "zoom",
        group.by     = "FinalCellType",
        molecules    = gene,
        nmols        = 20000,
        mols.size    = 0.3,
        mols.cols    = unname(marker_colors[gene]),
        border.color = "white",
        border.size  = 0.1,
        cols         = celltype_colors,
        alpha        = 0.3,
        coord.fixed  = FALSE,
        axes         = TRUE
      ) + ggtitle(paste(region_name, "-", gene)) +
        theme(legend.position = "right"),
      error = function(e) { message("  Plot failed for ", gene, ": ", e$message); NULL }
    )
    if (!is.null(p)) print(p)
  }

  dev.off()
  message("  Saved: ", out_pdf)

  region_obj
}

# ── Process each region ────────────────────────────────────────────────────────

# section 1: Crop bounds: x=[1437, 1802]  y=[1339, 1984]
ids_kr1 <- read_xenium_csv(file.path(REGION_DIR, "Xenium_KidneyRight_Selection_1_coordinates.csv"))
region_kr1 <- subset_region(ids_kr1, SAMPLE_RIGHT, IMAGE_RIGHT, "KidneyRight_Sel1")

# section 2: Crop bounds: x=[1922, 2475]  y=[1262, 1682]
ids_kr2 <- read_xenium_csv(file.path(REGION_DIR, "Xenium_KidneyRight_Selection_2_coordinates.csv"))
region_kr2 <- subset_region(ids_kr2, SAMPLE_RIGHT, IMAGE_RIGHT, "KidneyRight_Sel2")

# section 3: Crop bounds: x=[1244, 1665]  y=[1996, 2517]
ids_kr3 <- read_xenium_csv(file.path(REGION_DIR, "Xenium_KidneyRight_Selection_3_coordinates.csv"))
region_kr3 <- subset_region(ids_kr3, SAMPLE_RIGHT, IMAGE_RIGHT, "KidneyRight_Sel3")

# section 1: Crop bounds: x=[1775, 2408]  y=[1741, 2245]
ids_cp1 <- read_xenium_csv(file.path(REGION_DIR, "Xenium_KidneyLeft_Selection_1_coordinates.csv"))
region_cp1 <- subset_region(ids_cp1, SAMPLE_LEFT, IMAGE_LEFT, "KidneyLeft_Sel1")

# section 2: Crop bounds: x=[2277, 2859]  y=[1499, 1850]
ids_cp2 <- read_xenium_csv(file.path(REGION_DIR, "Xenium_KidneyLeft_Selection_2_coordinates.csv"))
region_cp2 <- subset_region(ids_cp2, SAMPLE_LEFT, IMAGE_LEFT, "KidneyLeft_Sel2")

# section 3: Crop bounds: x=[1457, 2122]  y=[2109, 2895]
ids_cp3 <- read_xenium_csv(file.path(REGION_DIR, "Xenium_KidneyLeft_Selection_3_coordinates.csv"))
region_cp3 <- subset_region(ids_cp3, SAMPLE_LEFT, IMAGE_LEFT, "KidneyLeft_Sel3")

message("\n==> Per-region processing complete.")

# ── Merge all non-NULL region objects ─────────────────────────────────────────
region_list <- list(
  KidneyRight_Sel1 = region_kr1,
  KidneyRight_Sel2 = region_kr2,
  KidneyRight_Sel3 = region_kr3,
  KidneyLeft_Sel1    = region_cp1,
  KidneyLeft_Sel2    = region_cp2,
  KidneyLeft_Sel3    = region_cp3
)
region_list <- Filter(Negate(is.null), region_list)

# ── Cell type fraction per region ───────────────────────────────────────────
# Uses the 6 individual per-region objects (before merge), since the merged
# object's "Region" column only distinguishes KidneyRight vs. KidneyLeft, not
# the 3 sub-selections within each.
celltype_fraction <- do.call(rbind, lapply(names(region_list), function(rn) {
  tab <- table(region_list[[rn]]$FinalCellType)
  data.frame(
    Region   = rn,
    CellType = names(tab),
    nCells   = as.integer(tab),
    Fraction = as.numeric(tab) / sum(tab)
  )
}))

write.csv(
  celltype_fraction,
  file.path(OUT_DIR, sprintf("CellTypeFraction_Xenium_SelectedRegions_%s.csv",
                              RequestedCategoriesName)),
  row.names = FALSE
)
message("  Saved cell type fraction table.")

celltype_fraction$CellType <- factor(celltype_fraction$CellType, levels = all_celltypes)
celltype_fraction$Region   <- factor(celltype_fraction$Region, levels = names(region_list))

FractionBarPlot <- ggplot(celltype_fraction, aes(x = Region, y = Fraction, fill = CellType)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = celltype_colors, drop = FALSE) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Fraction of cells", x = NULL, title = "Cell type composition by region")

ggsave(
  file.path(OUT_DIR, sprintf("CellTypeFraction_Xenium_SelectedRegions_%s.pdf",
                              RequestedCategoriesName)),
  FractionBarPlot, width = 6, height = 6
)
message("  Saved cell type fraction barplot.")

# ── Spatial proximity: Uro cells to macrophage/fibroblast/smooth muscle/PC ──
# Nearest-neighbor centroid distance (tissue microns) from each urothelial
# cell to the closest cell of each target type, per region. Uses the cropped
# "zoom" FOV centroids (same coordinate space used for the ImageDimPlot
# figures). "Smooth muscle" maps to "Peri/VSMC" — the manual annotation never
# separated pericytes from vascular smooth muscle into their own clusters.
proximity_targets <- c(Macrophage = "Mac", Fibroblast = "Fib",
                        SmoothMuscle = "Peri/VSMC", TCell = "T lymph", PrincipalCell = "PC", EndothelialCell= "Endo")

# Brute-force min distance (fine at region-crop scale, hundreds-few thousand
# cells); returns NA for a source cell when the target type isn't present.
nearest_dist <- function(from_xy, to_xy) {
  if (nrow(from_xy) == 0) return(numeric(0))
  if (nrow(to_xy) == 0) return(rep(NA_real_, nrow(from_xy)))
  d <- sqrt(outer(from_xy[, 1], to_xy[, 1], "-")^2 + outer(from_xy[, 2], to_xy[, 2], "-")^2)
  apply(d, 1, min)
}

proximity_df <- do.call(rbind, lapply(names(region_list), function(rn) {
  robj <- region_list[[rn]]

  coords <- GetTissueCoordinates(robj, image = "zoom", which = "centroids")
  rownames(coords) <- coords$cell
  coords <- coords[Cells(robj), ]

  celltype <- as.character(robj$FinalCellType)
  uro_idx  <- which(celltype == "Uro")
  if (length(uro_idx) == 0) {
    message(sprintf("  [skip] %s: no Uro cells", rn))
    return(NULL)
  }
  uro_xy    <- as.matrix(coords[uro_idx, c("x", "y")])
  uro_cells <- Cells(robj)[uro_idx]

  do.call(rbind, lapply(names(proximity_targets), function(target_name) {
    target_type <- proximity_targets[[target_name]]
    target_idx  <- which(celltype == target_type)
    target_xy   <- as.matrix(coords[target_idx, c("x", "y")])
    data.frame(
      Region     = rn,
      UroCell    = uro_cells,
      TargetType = target_name,
      nTarget    = length(target_idx),
      Distance   = nearest_dist(uro_xy, target_xy)
    )
  }))
}))

write.csv(
  proximity_df,
  file.path(OUT_DIR, sprintf("UroProximity_Xenium_SelectedRegions_%s.csv",
                              RequestedCategoriesName)),
  row.names = FALSE
)
message("  Saved Uro proximity distance table.")

proximity_summary <- proximity_df %>%
  group_by(Region, TargetType) %>%
  summarise(
    nUro       = sum(!is.na(Distance)),
    nTarget    = dplyr::first(nTarget),
    MedianDist = median(Distance, na.rm = TRUE),
    MeanDist   = mean(Distance, na.rm = TRUE),
    .groups    = "drop"
  )
write.csv(
  proximity_summary,
  file.path(OUT_DIR, sprintf("UroProximitySummary_Xenium_SelectedRegions_%s.csv",
                              RequestedCategoriesName)),
  row.names = FALSE
)
message("  Saved Uro proximity summary table.")

proximity_df$Region <- factor(proximity_df$Region, levels = names(region_list))
ProximityBoxplot <- ggplot(proximity_df, aes(x = TargetType, y = Distance, fill = TargetType)) +
  geom_boxplot(outlier.size = 0.5, na.rm = TRUE) +
  facet_wrap(~ Region, nrow = 2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(y = "Nearest-neighbor distance from Uro cell (µm)", x = NULL,
       title = "Distance from urothelial cells to nearby cell types")

ggsave(
  file.path(OUT_DIR, sprintf("UroProximityBoxplot_Xenium_SelectedRegions_%s.pdf",
                              RequestedCategoriesName)),
  ProximityBoxplot, width = 6, height = 6
)
message("  Saved Uro proximity boxplot.")

message(sprintf("\n==> Merging %d non-empty region(s) ...", length(region_list)))

if (length(region_list) < 2) {
  message("[warning] Fewer than 2 regions available — skipping merge and downstream analysis.")
  quit(save = "no", status = 0)
}

UrotheliumSelectedRegions <- merge(
  region_list[[1]],
  y            = region_list[-1],
  add.cell.ids = names(region_list),
  project      = "XeniumUrotheliumSelectedRegions"
)
# Region label: first two components of the cell name after merge (e.g.
# "KidneyLeft_Sel3", not just "KidneyLeft") — add.cell.ids prefixes every cell
# with "<region_list name>_", and region_list names are all "Kidney{Left,
# Right}_Sel{1,2,3}", i.e. two underscore-separated parts.
UrotheliumSelectedRegions$Region <- sapply(
  strsplit(colnames(UrotheliumSelectedRegions), "_"),
  function(parts) paste(parts[1:2], collapse = "_")
)

names(region_list)

UrotheliumSelectedRegions$Region <- factor(UrotheliumSelectedRegions$Region, levels = names(region_list))

Requested_MARKERS <- c(
  "Aqp2", "Scnn1g", # principal cells
  "Col1a2", "Col5a2", # Fibroblast
  "Acta2", "Myh11", # VSMC
  "Emcn", "Kdr", # Endo
  "Csf1r","Il34", # Mac
  "Il2", "Il2ra" # T cell
)
# generate the dotplot for these nearby cell markers for these six regions
  DotPlotUro <- DotPlot(
    UrotheliumSelectedRegions,
    features  = Requested_MARKERS,
    group.by  = "Region",
    cols      = c("lightgrey", "red"),
    dot.scale = 8
  ) + RotatedAxis() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    guides(
      color = guide_colorbar(title = "AveExp"),
      size  = guide_legend(title  = "PerExp")
    )

  ggsave(
    file.path(OUT_DIR, sprintf("DotPlot_Xenium_UrotheliumSelectedRegions_%s.pdf",
                                RequestedCategoriesName)),
    DotPlotUro, width = 6, height = 3.5
  )
  message("  Saved DotPlot.")


out_merged <- file.path(OUT_DIR, "Xenium_UrotheliumSelectedSixRegions.rds")
saveRDS(UrotheliumSelectedRegions, out_merged)
message("  Saved merged regions object: ", out_merged)
UrotheliumSelectedRegions <-readRDS(out_merged)
# ── Subset to urothelium cells (RCTD dominant cell type) ──────────────────────
# For Xenium single cells, use rctd_dominant_celltype rather than cluster number.
# If you prefer cluster-based: subset(UrotheliumSelectedRegions, idents = "9")
OnlyUrotheliumCells <- subset(UrotheliumSelectedRegions,
                               subset = FinalCellType == "Uro")
message(sprintf("  Urothelium cells in selected regions: %d", ncol(OnlyUrotheliumCells)))

OnlyUrotheliumCells%>% head()

# ── DotPlot: marker expression by region ──────────────────────────────────────
markers_in_uro <- c( "Upk1b", "Krt5", "Krt14","Krt15" ,"Krt19","Krt20","Cd24a","Spp1", "Havcr1",  "Vcam1", "Cd74","Fosl1" )

if (length(markers_in_uro) > 0 && ncol(OnlyUrotheliumCells) > 0) {
  DotPlotUro <- DotPlot(
    OnlyUrotheliumCells,
    features  = markers_in_uro,
    group.by  = "Region",
    cols      = c("lightgrey", "red"),
    dot.scale = 8
  ) + RotatedAxis() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    guides(
      color = guide_colorbar(title = "AveExp"),
      size  = guide_legend(title  = "PerExp")
    )

  ggsave(
    file.path(OUT_DIR, sprintf("DotPlot_Xenium_UrotheliumSelectedRegions_%s.pdf",
                                RequestedCategoriesName)),
    DotPlotUro, width = 6, height = 3.5
  )
  message("  Saved DotPlot.")
}

# ── DEG analysis: region-specific markers within urothelium cells ──────────────
OnlyUrotheliumCells <- JoinLayers(OnlyUrotheliumCells)
Idents(OnlyUrotheliumCells) <- OnlyUrotheliumCells$Region

RegionalDEGs <- FindAllMarkers(
  OnlyUrotheliumCells,
  only.positive   = TRUE,
  min.pct         = 0.2,
  logfc.threshold = 0.2,
  test.use        = "wilcox"
)

write.csv(
  RegionalDEGs,
  file.path(OUT_DIR, sprintf("DEGs_Xenium_UrotheliumSelectedRegions_%s.csv",
                              RequestedCategoriesName)),
  row.names = FALSE
)
message("  Saved DEG table.")

Top5MarkersPerRegion <- RegionalDEGs %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

if (nrow(Top5MarkersPerRegion) > 0) {
  Top5Heatmap <- DoHeatmap(
    OnlyUrotheliumCells,
    features = Top5MarkersPerRegion$gene,
    group.by = "Region",
    size     = 3
  ) + scale_fill_gradientn(
    colors   = c("blue", "white", "red"),
    na.value = "white"
  ) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(
    file.path(OUT_DIR, sprintf("Heatmap_Top5DEGs_Xenium_UrotheliumSelectedRegions_%s.pdf",
                                RequestedCategoriesName)),
    Top5Heatmap, width = 6, height = 4
  )
  message("  Saved DEG heatmap.")
}

message("\n==> All done.")


# We only selected these two healthy Xenium slides and check the Biomarkers
object@meta.data %>% head()
HealthyObject <- subset(object, sample_id %in% c("ShamL","ShamR"))

# Generate the DotPlot for the cell type markers, each please choose two 
CellTypeMarkers <- c(
  "Slc5a2", "Slc5a12",      # PT-S1
  "Slc22a6", "Slc13a3",     # PT-S2/S3
  "Cyp7b1",      # PT-S3
  "Slc7a13",      # PT-S3/TAL
  "Havcr1",      # Inj-PT
  "Vcam1",          # FR-PT
  "Mki67",      # Prolif-Inj-PT
  "Slc14a2", "Aqp1",        # DTL
  "Slc12a1", "Umod",        # TAL
  "Slc12a3", "Trpm6",       # DCT
  "Aqp2", "Scnn1g",           # PC
  "Slc4a1", "Aqp6",         # A-IC
  "Slc26a4", "Slc4a9",      # B-IC
  "Nphs2", "Synpo",        # Podo
  "Ehd3",          # Podo/GEnC
  "Upk1b", "Krt5",          # Uro
  "Kdr", "Emcn",            # Endo
  "Pde3a", "Acta2", "Myh11",         # Peri/VSMC
  "Col1a1", "Col1a2",        # Fib
  "Csf1r", "Cd74" ,          # Mac
  "Il2","Il2ra"
)

Idents(HealthyObject) <- factor(HealthyObject$FinalCellType, levels = rev(all_celltypes))
 
p <- DotPlot(
  HealthyObject,
  features = CellTypeMarkers,
  cols     = c("lightgrey", "blue"),
  dot.scale = 6
) +
  RotatedAxis() +
    guides(
      color = guide_colorbar(title = "AveExp"),
      size  = guide_legend(title  = "PerExp")
    )
 
print(p)
 
## Save
ggsave(file.path(OUT_DIR,"Xenium_dotplot_celltypes.pdf"), plot = p, width = 10, height = 8)


# we only selected the Urothelium cells

HealthyUrotheliumObject <- subset(HealthyObject, FinalCellType == "Uro")

Idents(HealthyUrotheliumObject)<- "sample_id"
# check the expression of 
markers_in_uro <- c( "Upk1b", "Krt5", "Krt14","Krt15" ,"Krt19","Krt20","Cd24a","Spp1", "Havcr1",  "Vcam1", "Cd74","Fosl1" )
DotPlot(HealthyUrotheliumObject, features = markers_in_uro)

p <- DotPlot(
  HealthyUrotheliumObject,
  features = markers_in_uro,
  cols     = c("lightgrey", "blue"),
  #dot.scale = 6
) +
  RotatedAxis() +  coord_flip()+
    guides(
      color = guide_colorbar(title = "AveExp"),
      size  = guide_legend(title  = "PerExp")
    ) 
 
print(p)
 
## Save
ggsave(file.path(OUT_DIR,"Xenium_dotplot_uroinjurymarker_onlyinurothelium.pdf"), plot = p, width = 3, height = 5)
