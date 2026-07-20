################################################################################
# 06_VisiumHD_Niche9_NeighborhoodAnalysis.R
#
# Spatial neighborhood cell-type composition around urothelium bins in the
# VisiumHD 8 um-bin data — niche 9 (seurat_cluster.harmony.projected == "9")
# bins further restricted to rctd_dominant_celltype == "Urothelium", so that
# non-urothelial bins that merely fall in niche 9 don't dilute the "urothelium
# neighborhood" signal.
#
# For every niche-9 bin, all OTHER bins on the same image (kidney3p or
# kidney5p — coordinates are never mixed across images) within a fixed
# physical radius are treated as spatial neighbors. Each niche-9 bin's
# "neighborhood composition" is the mean RCTD weight vector of its
# neighbors — i.e. what cell types physically surround that urothelium bin,
# not the bin's own deconvolution call.
#
# Radius: 50 um (~6 bin-widths; approximates the immediate subepithelial
# microenvironment). Converted from physical um to the image's pixel space
# using each VisiumV2 image's own scale.factors$spot (pixels per 8 um bin
# edge), so the two kidney images (which can have different scale factors)
# are each handled in their native pixel scale.
#
# Same RCTD cell-type collapsing as 05_VisiumHD_Niche9_Deconvolution_Matrix.R:
#   - Asc-Vasa-Recta / Desc-Vasa-Recta / Vas-Efferens / Vas-Afferens summed
#     into "Endo".
#   - "Per" relabeled "Peri/VSMC".
#
# Two groupings, mirroring 05_VisiumHD_Niche9_Deconvolution_Matrix.R:
#   (A) all niche-9 bins (both kidneys), grouped by sample_id
#   (B) bins within the five manually-selected urothelium regions from
#       04_VisiumHD_Urothelium_RegionSubsets.R that are themselves annotated
#       rctd_dominant_celltype == "Urothelium" (independent of niche 9 cluster
#       membership), grouped by Region
#
# Input:  VisiumHD_kidney_both_deconvolved.rds
#         Kidney1TopUroSection.csv / Kidney1MiddleUroSection.csv / Kidney1BottomUroSection.csv
#         Kidney3P_LeftUrothelium.csv / Kidney3P_RightUrothelium.csv
# Output: NeighborhoodMatrix_Niche9_AllBins.csv / NeighborhoodMatrix_Niche9_MeanBySample.csv
#         Heatmap_Niche9_Neighborhood_MeanBySample.pdf / Barplot_Niche9_Neighborhood_MeanBySample.pdf
#         NeighborhoodMatrix_Niche9_FiveRegions.csv / NeighborhoodMatrix_Niche9_MeanByRegion.csv
#         Heatmap_Niche9_Neighborhood_MeanByRegion.pdf / Barplot_Niche9_Neighborhood_MeanByRegion.pdf
#         UroProximity_Niche9_DistanceToCellTypes.csv / UroProximitySummary_Niche9_DistanceToCellTypes.csv
#         UroProximityBoxplot_Niche9_DistanceToCellTypes.pdf
#         UroProximity_Niche9_DistanceToFixedCellTypes.csv / UroProximitySummary_Niche9_DistanceToFixedCellTypes.csv
#         UroProximityBoxplot_Niche9_DistanceToFixedCellTypes.pdf
#         UroProximity_Niche9_DistanceToFixedCellTypes_ByRegion.csv
#         UroProximitySummary_Niche9_DistanceToFixedCellTypes_ByRegion.csv
#         UroProximityBoxplot_Niche9_DistanceToFixedCellTypes_ByRegion.pdf
################################################################################

library(Seurat)
library(RANN)
library(ggplot2)
library(dplyr)
library(tidyr)

# ── Paths ─────────────────────────────────────────────────────────────────────
HD_BASE <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSpatialData/VisiumHD"
OUT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/SpatialScripts/output"

OBJECT_RDS <- file.path(OUT_DIR, "VisiumHD_kidney_both_deconvolved.rds")
NICHE_ID          <- "9"
NEIGHBOR_RADIUS_UM <- 50
BIN_WIDTH_UM       <- 8
KNN_K              <- 300  # generous upper bound on bins within NEIGHBOR_RADIUS_UM (~122 expected)

# ── Collapse related RCTD cell types for interpretability ─────────────────────
VASCULAR_SUBTYPES <- c("rctd_Asc-Vasa-Recta", "rctd_Desc-Vasa-Recta",
                        "rctd_Vas-Efferens", "rctd_Vas-Afferens")

collapse_rctd_celltypes <- function(df) {
  vascular_present <- intersect(VASCULAR_SUBTYPES, colnames(df))
  if (length(vascular_present) > 0) {
    df$rctd_Endo <- rowSums(df[, c("rctd_Endo", vascular_present), drop = FALSE])
    df <- df[, !(colnames(df) %in% vascular_present)]
  }
  if ("rctd_Per" %in% colnames(df)) {
    colnames(df)[colnames(df) == "rctd_Per"] <- "rctd_Peri/VSMC"
  }
  if ("rctd_dominant_celltype" %in% colnames(df)) {
    vascular_labels <- sub("^rctd_", "", VASCULAR_SUBTYPES)
    df$rctd_dominant_celltype[df$rctd_dominant_celltype %in% vascular_labels] <- "Endo"
    df$rctd_dominant_celltype[df$rctd_dominant_celltype == "Per"] <- "Peri/VSMC"
  }
  df
}

# ── Load deconvolved object ────────────────────────────────────────────────────
message("==> Loading VisiumHD deconvolved object ...")
object <- readRDS(OBJECT_RDS)

meta <- object@meta.data
meta <- collapse_rctd_celltypes(meta)
rctd_cols <- grep("^rctd_", colnames(meta), value = TRUE)
rctd_cols <- setdiff(rctd_cols, "rctd_dominant_celltype")
message(sprintf("RCTD cell types after collapsing: %d", length(rctd_cols)))
message("  ", paste(sub("^rctd_", "", rctd_cols), collapse = ", "))

# ── Helper: mean neighbor-composition matrix for a set of target bin IDs ──────
# For each target bin, average the (collapsed) rctd_* weights of every OTHER
# bin on the same image within NEIGHBOR_RADIUS_UM. Images are never mixed.
compute_neighborhood_composition <- function(object, meta, rctd_cols, target_ids,
                                              radius_um = NEIGHBOR_RADIUS_UM) {
  all_images   <- Images(object)
  result_list  <- list()

  for (img in all_images) {
    coords <- GetTissueCoordinates(object, image = img, which = "centroids")
    rownames(coords) <- coords$cell
    img_ids      <- rownames(coords)
    img_targets  <- intersect(target_ids, img_ids)
    if (length(img_targets) == 0) next

    spot_px   <- object[[img]]@scale.factors$spot  # pixels per BIN_WIDTH_UM bin edge
    radius_px <- (radius_um / BIN_WIDTH_UM) * spot_px
    message(sprintf("  Image %s: spot=%.2f px/bin, radius=%.1f px (%d target bins)",
                    img, spot_px, radius_px, length(img_targets)))

    coord_mat <- as.matrix(coords[, c("x", "y")])
    query_mat <- coord_mat[img_targets, , drop = FALSE]
    k <- min(KNN_K, nrow(coord_mat))

    nn <- RANN::nn2(coord_mat, query_mat, k = k)
    if (any(nn$nn.dists[, k] <= radius_px)) {
      warning(sprintf(
        "  KNN_K=%d may be too small for image %s — some bins have >%d ",
        k, img, k), "neighbors within the radius; true neighbor counts may be truncated.")
    }

    weight_mat <- as.matrix(meta[img_ids, rctd_cols])

    comp <- matrix(NA_real_, nrow = length(img_targets), ncol = length(rctd_cols),
                   dimnames = list(img_targets, rctd_cols))
    n_neighbors <- integer(length(img_targets))

    for (i in seq_along(img_targets)) {
      nbr_ids  <- img_ids[nn$nn.idx[i, ]]
      nbr_dist <- nn$nn.dists[i, ]
      keep <- (nbr_dist <= radius_px) & (nbr_ids != img_targets[i])
      n_neighbors[i] <- sum(keep)
      if (n_neighbors[i] > 0) {
        comp[i, ] <- colMeans(weight_mat[nbr_ids[keep], , drop = FALSE], na.rm = TRUE)
      }
    }

    comp_df <- as.data.frame(comp)
    comp_df$n_neighbors <- n_neighbors
    result_list[[img]] <- comp_df
  }

  do.call(rbind, result_list)
}

# ── Helper: heatmap + barplot from a (group x celltype) mean weight matrix ────
plot_composition <- function(mean_mat, group_col_name, out_prefix) {
  mean_df <- as.data.frame(mean_mat)
  mean_df[[group_col_name]] <- rownames(mean_df)
  long_df <- pivot_longer(mean_df, -all_of(group_col_name),
                           names_to = "celltype", values_to = "weight")
  long_df$celltype <- sub("^rctd_", "", long_df$celltype)

  n_groups <- length(unique(long_df[[group_col_name]]))

  p_heat <- ggplot(long_df, aes(x = .data[[group_col_name]], y = celltype, fill = weight)) +
    geom_tile(color = "white") +
    scale_fill_gradientn(colors = c("white", "#FED98E", "#FE8929", "#CC4C02", "#7F0000"),
                          name = "Mean neighbor\nRCTD weight") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = NULL)
  ggsave(file.path(OUT_DIR, sprintf("Heatmap_%s.pdf", out_prefix)), p_heat,
         width = max(5, 1 + n_groups * 1.2), height = 8)

  p_bar <- ggplot(long_df, aes(x = .data[[group_col_name]], y = weight, fill = celltype)) +
    geom_col(position = "fill") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = "Neighborhood cell type proportion", fill = "Cell type")
  ggsave(file.path(OUT_DIR, sprintf("Barplot_%s.pdf", out_prefix)), p_bar,
         width = max(5, 2 + n_groups * 1.2), height = 6)

  message("  Saved: Heatmap_", out_prefix, ".pdf and Barplot_", out_prefix, ".pdf")
}

################################################################################
# Part A: niche-9 neighborhoods across the whole object, grouped by sample_id
################################################################################
message("\n==> Part A: niche 9 neighborhood composition, grouped by sample_id ...")

niche9_ids <- rownames(meta)[
  meta$seurat_cluster.harmony.projected == NICHE_ID &
    meta$rctd_dominant_celltype %in% "Urothelium"
]
message(sprintf("  Niche %s bins with rctd_dominant_celltype == Urothelium: %d", NICHE_ID, length(niche9_ids)))

niche9_neighborhood <- compute_neighborhood_composition(object, meta, rctd_cols, niche9_ids)
niche9_neighborhood$sample_id <- meta[rownames(niche9_neighborhood), "sample_id"]

write.csv(niche9_neighborhood, file.path(OUT_DIR, "NeighborhoodMatrix_Niche9_AllBins.csv"))

mean_by_sample <- niche9_neighborhood %>%
  group_by(sample_id) %>%
  summarise(across(all_of(rctd_cols), ~ mean(.x, na.rm = TRUE))) %>%
  as.data.frame()
rownames(mean_by_sample) <- mean_by_sample$sample_id
mean_by_sample$sample_id <- NULL
write.csv(mean_by_sample, file.path(OUT_DIR, "NeighborhoodMatrix_Niche9_MeanBySample.csv"))

plot_composition(mean_by_sample, "sample_id", "Niche9_Neighborhood_MeanBySample")

################################################################################
# Part B: niche-9 neighborhoods within the five integrated urothelium regions
################################################################################
message("\n==> Part B: niche 9 neighborhood composition within the five integrated regions ...")

region_specs <- list(
  Pelviscalyceal_urothelium   = list(csv = "Kidney1TopUroSection.csv",     prefix = "kidney_"),
  PelvisCalyxLumen_urothelium = list(csv = "Kidney1MiddleUroSection.csv",  prefix = "kidney_"),
  Periurothelial_urothelium   = list(csv = "Kidney1BottomUroSection.csv",  prefix = "kidney_"),
  ThreeSideLeftUrothelium     = list(csv = "Kidney3P_LeftUrothelium.csv",  prefix = "kidney3p_"),
  ThreeSideRightUrothelium    = list(csv = "Kidney3P_RightUrothelium.csv", prefix = "kidney3p_")
)

region_id_list <- list()
for (region_name in names(region_specs)) {
  spec     <- region_specs[[region_name]]
  csv      <- read.csv(file.path(HD_BASE, spec$csv))
  cell_ids <- paste0(spec$prefix, csv$Barcode)
  matched  <- intersect(cell_ids, rownames(meta))
  message(sprintf("  %-28s requested: %5d | matched: %5d",
                  region_name, length(cell_ids), length(matched)))
  if (length(matched) == 0) {
    stop(sprintf(
      "No cells matched for region '%s' (prefix '%s'). ",
      region_name, spec$prefix),
      "Barcode prefix or CSV is likely wrong — do not proceed silently."
    )
  }
  region_id_list[[region_name]] <- matched
}

# Restrict each region's bins to those annotated as Urothelium
# (rctd_dominant_celltype == "Urothelium"), independent of niche 9 cluster
# membership — some urothelium-annotated bins fall outside cluster 9, and
# those are included here (unlike niche9_ids used in Part A).
region_uro_ids <- lapply(region_id_list, function(ids) {
  ids[meta[ids, "rctd_dominant_celltype"] %in% "Urothelium"]
})
message("  Urothelium-annotated bins per region:")
print(sapply(region_uro_ids, length))

region_uro_ids_all <- unlist(region_uro_ids, use.names = FALSE)
region_uro_neighborhood <- compute_neighborhood_composition(object, meta, rctd_cols, region_uro_ids_all)

five_region_neighborhood_list <- lapply(names(region_uro_ids), function(region_name) {
  ids <- region_uro_ids[[region_name]]
  df  <- region_uro_neighborhood[ids, , drop = FALSE]
  df$Region <- region_name
  df
})
five_region_neighborhood <- do.call(rbind, five_region_neighborhood_list)

write.csv(five_region_neighborhood, file.path(OUT_DIR, "NeighborhoodMatrix_Niche9_FiveRegions.csv"))

mean_by_region <- five_region_neighborhood %>%
  group_by(Region) %>%
  summarise(across(all_of(rctd_cols), ~ mean(.x, na.rm = TRUE))) %>%
  as.data.frame()
rownames(mean_by_region) <- mean_by_region$Region
mean_by_region$Region <- NULL
write.csv(mean_by_region, file.path(OUT_DIR, "NeighborhoodMatrix_Niche9_MeanByRegion.csv"))

plot_composition(mean_by_region, "Region", "Niche9_Neighborhood_MeanByRegion")

################################################################################
# Part C: nearest-neighbor distance from Uro bins (niche 9,
# rctd_dominant_celltype == "Urothelium") to the closest bin whose own
# rctd_dominant_celltype is a given target type
################################################################################
message("\n==> Part C: Uro bin proximity to nearby cell types ...")

# Originally-requested target types, kept as a fixed baseline set.
fixed_targets <- c(Endothelial = "Endo", Fibroblast = "Fib", Macrophage = "Macro",
                    SmoothMuscle = "Peri/VSMC", PrincipalCell = "PC", TCell = "T lymph",CDtrans= "CD-Trans")

# All non-Urothelium rctd_dominant_celltype calls present in the object are
# candidate targets; distance is computed to every one of them below, then
# narrowed to the 5 with the smallest median distance from Uro bins — i.e.
# the 5 cell types that are actually spatially nearest — and combined with
# the fixed baseline set above.
candidate_target_types <- setdiff(unique(meta$rctd_dominant_celltype), c(NA, "Urothelium"))
candidate_targets <- setNames(candidate_target_types, candidate_target_types)

# Distance is converted from pixels to physical um via each image's own
# scale.factors$spot (pixels per BIN_WIDTH_UM bin edge), same conversion used
# for NEIGHBOR_RADIUS_UM above. Images are never mixed.
compute_uro_proximity <- function(object, meta, uro_ids, targets) {
  all_images  <- Images(object)
  result_list <- list()

  for (img in all_images) {
    coords <- GetTissueCoordinates(object, image = img, which = "centroids")
    rownames(coords) <- coords$cell
    img_ids <- rownames(coords)
    img_uro <- intersect(uro_ids, img_ids)
    if (length(img_uro) == 0) next

    spot_px   <- object[[img]]@scale.factors$spot
    px_to_um  <- BIN_WIDTH_UM / spot_px
    coord_mat <- as.matrix(coords[, c("x", "y")])
    uro_mat   <- coord_mat[img_uro, , drop = FALSE]

    for (target_name in names(targets)) {
      target_type <- targets[[target_name]]
      target_ids  <- img_ids[meta[img_ids, "rctd_dominant_celltype"] %in% target_type]
      target_ids  <- setdiff(target_ids, img_uro)  # a bin can't be its own nearest neighbor
      if (length(target_ids) == 0) {
        message(sprintf("  [skip] image %s, target %s: no bins with dominant type '%s'",
                        img, target_name, target_type))
        dist_um <- rep(NA_real_, length(img_uro))
      } else {
        target_mat <- coord_mat[target_ids, , drop = FALSE]
        nn <- RANN::nn2(target_mat, uro_mat, k = 1)
        dist_um <- nn$nn.dists[, 1] * px_to_um
      }
      result_list[[paste(img, target_name)]] <- data.frame(
        Image      = img,
        UroBin     = img_uro,
        TargetType = target_name,
        nTarget    = length(target_ids),
        Distance   = dist_um
      )
    }
  }

  do.call(rbind, result_list)
}

uro_proximity_all <- compute_uro_proximity(object, meta, niche9_ids, candidate_targets)
uro_proximity_all$sample_id <- meta[uro_proximity_all$UroBin, "sample_id"]

# Rank candidate target types by overall median distance from Uro bins and
# keep the 5 nearest.
nearest_ranking <- uro_proximity_all %>%
  group_by(TargetType) %>%
  summarise(MedianDist = median(Distance, na.rm = TRUE), .groups = "drop") %>%
  arrange(MedianDist)
top10_nearest <- head(nearest_ranking$TargetType, 10)
message("  Top 10 nearest target cell types: ", paste(top10_nearest, collapse = ", "))

# Combine the fixed baseline set with any top-10-nearest types not already in
# it, keeping the friendly display names for the fixed set and using the raw
# rctd_dominant_celltype label as the display name for the rest.
extra_nearest   <- setdiff(top10_nearest, fixed_targets)
combined_targets <- c(fixed_targets, setNames(extra_nearest, extra_nearest))
message("  Combined target set: ", paste(names(combined_targets), collapse = ", "))

display_name_by_type <- setNames(names(combined_targets), combined_targets)
uro_proximity <- uro_proximity_all %>%
  filter(TargetType %in% combined_targets) %>%
  mutate(TargetType = unname(display_name_by_type[TargetType]))

write.csv(uro_proximity, file.path(OUT_DIR, "UroProximity_Niche9_DistanceToCellTypes.csv"),
          row.names = FALSE)
message("  Saved Uro proximity distance table.")

uro_proximity_summary <- uro_proximity %>%
  group_by(sample_id, TargetType) %>%
  summarise(
    nUro       = sum(!is.na(Distance)),
    nTarget    = dplyr::first(nTarget),
    MedianDist = median(Distance, na.rm = TRUE),
    MeanDist   = mean(Distance, na.rm = TRUE),
    .groups    = "drop"
  )
write.csv(uro_proximity_summary,
          file.path(OUT_DIR, "UroProximitySummary_Niche9_DistanceToCellTypes.csv"),
          row.names = FALSE)
message("  Saved Uro proximity summary table.")

uro_proximity$TargetType <- factor(uro_proximity$TargetType, levels = names(combined_targets))
UroProximityBoxplot <- ggplot(uro_proximity, aes(x = TargetType, y = Distance, fill = TargetType)) +
  geom_boxplot(outlier.size = 0.5, na.rm = TRUE) +
  facet_wrap(~ sample_id, ncol = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(y = "Nearest-neighbor distance from Uro bin (µm)", x = NULL)

ggsave(file.path(OUT_DIR, "UroProximityBoxplot_Niche9_DistanceToCellTypes.pdf"),
       UroProximityBoxplot, width = 3, height = 4)
message("  Saved Uro proximity boxplot.")

################################################################################
# Part D: same proximity analysis restricted to fixed_targets only (no
# top-nearest ranking/combination) — a separate, simpler distance plot.
################################################################################
message("\n==> Part D: Uro bin proximity to fixed_targets only ...")

uro_proximity_fixed <- compute_uro_proximity(object, meta, niche9_ids, fixed_targets)
uro_proximity_fixed$sample_id <- meta[uro_proximity_fixed$UroBin, "sample_id"]

write.csv(uro_proximity_fixed,
          file.path(OUT_DIR, "UroProximity_Niche9_DistanceToFixedCellTypes.csv"),
          row.names = FALSE)
message("  Saved Uro proximity (fixed_targets) distance table.")

uro_proximity_fixed_summary <- uro_proximity_fixed %>%
  group_by(sample_id, TargetType) %>%
  summarise(
    nUro       = sum(!is.na(Distance)),
    nTarget    = dplyr::first(nTarget),
    MedianDist = median(Distance, na.rm = TRUE),
    MeanDist   = mean(Distance, na.rm = TRUE),
    .groups    = "drop"
  )
write.csv(uro_proximity_fixed_summary,
          file.path(OUT_DIR, "UroProximitySummary_Niche9_DistanceToFixedCellTypes.csv"),
          row.names = FALSE)
message("  Saved Uro proximity (fixed_targets) summary table.")

uro_proximity_fixed$TargetType <- factor(uro_proximity_fixed$TargetType, levels = names(fixed_targets))
UroProximityBoxplotFixed <- ggplot(uro_proximity_fixed, aes(x = TargetType, y = Distance, fill = TargetType)) +
  geom_boxplot(outlier.size = 0.5, na.rm = TRUE) +
  facet_wrap(~ sample_id) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(y = "Nearest-neighbor distance from Uro bin (µm)", x = NULL)

ggsave(file.path(OUT_DIR, "UroProximityBoxplot_Niche9_DistanceToFixedCellTypes.pdf"),
       UroProximityBoxplotFixed, width = 5, height = 3)
message("  Saved Uro proximity (fixed_targets) boxplot.")

################################################################################
# Part E: same fixed_targets proximity analysis as Part D, but for the
# Urothelium-annotated bins within the five integrated regions (region_uro_ids
# from Part B), broken down by Region instead of sample_id.
################################################################################
message("\n==> Part E: Uro bin proximity to fixed_targets, by region ...")

uro_proximity_fixed_region <- compute_uro_proximity(object, meta, region_uro_ids_all, fixed_targets)

region_lookup <- setNames(rep(names(region_uro_ids), lengths(region_uro_ids)),
                           unlist(region_uro_ids, use.names = FALSE))
uro_proximity_fixed_region$Region <- region_lookup[uro_proximity_fixed_region$UroBin]

write.csv(uro_proximity_fixed_region,
          file.path(OUT_DIR, "UroProximity_Niche9_DistanceToFixedCellTypes_ByRegion.csv"),
          row.names = FALSE)
message("  Saved Uro proximity (fixed_targets, by region) distance table.")

uro_proximity_fixed_region_summary <- uro_proximity_fixed_region %>%
  group_by(Region, TargetType) %>%
  summarise(
    nUro       = sum(!is.na(Distance)),
    nTarget    = dplyr::first(nTarget),
    MedianDist = median(Distance, na.rm = TRUE),
    MeanDist   = mean(Distance, na.rm = TRUE),
    .groups    = "drop"
  )
write.csv(uro_proximity_fixed_region_summary,
          file.path(OUT_DIR, "UroProximitySummary_Niche9_DistanceToFixedCellTypes_ByRegion.csv"),
          row.names = FALSE)
message("  Saved Uro proximity (fixed_targets, by region) summary table.")

uro_proximity_fixed_region$TargetType <- factor(uro_proximity_fixed_region$TargetType, levels = names(fixed_targets))
uro_proximity_fixed_region$Region <- factor(uro_proximity_fixed_region$Region, levels = names(region_uro_ids))
UroProximityBoxplotFixedByRegion <- ggplot(uro_proximity_fixed_region,
                                            aes(x = TargetType, y = Distance, fill = Region)) +
  geom_boxplot(outlier.size = 0.2, na.rm = TRUE) +
  # facet_wrap(~ Region, ncol = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right") +
  labs(y = "Nearest-neighbor distance from Uro bin (µm)", x = NULL)

ggsave(file.path(OUT_DIR, "UroProximityBoxplot_Niche9_DistanceToFixedCellTypes_ByRegion.pdf"),
       UroProximityBoxplotFixedByRegion, width = 6, height = 3)
message("  Saved Uro proximity (fixed_targets, by region) boxplot.")

message("\n==> Done.")
