################################################################################
# 05_VisiumHD_Niche9_Deconvolution_Matrix.R
#
# Extract RCTD cell-type deconvolution weights (rctd_*) for niche 9 (the
# urothelium niche, seurat_cluster.harmony.projected == "9") and summarize
# them as a composition matrix, heatmap, and cell-type barplot — for:
#   (A) niche 9 across the whole deconvolved object (both kidneys, grouped by
#       sample_id: Visium_HD_3prime_Kidney vs Visium_HD_Kidney)
#   (B) niche 9 restricted to the five manually-selected urothelium regions
#       from 04_VisiumHD_Urothelium_RegionSubsets.R (OnlyUrotheliumCluster),
#       grouped by Region
#
# Before building either matrix, related RCTD cell types are collapsed:
#   - Asc-Vasa-Recta, Desc-Vasa-Recta, Vas-Efferens, Vas-Afferens weights are
#     summed into the existing "Endo" column (minor vasculature -> general
#     endothelium), then dropped.
#   - "Per" is relabeled "Peri/VSMC" (RCTD's reference calls both pericytes
#     and vascular smooth muscle cells "Per").
#
# Region barcode CSVs/prefixes mirror 04_VisiumHD_Urothelium_RegionSubsets.R.
# Only meta.data (not full subset Seurat objects) is pulled for part (B) since
# the rctd_* weights and cluster call are all that's needed — this avoids the
# heavy subset()/JoinLayers() calls script 04 uses for spatial marker plots.
#
# Input:  VisiumHD_kidney_both_deconvolved.rds
#         Kidney1TopUroSection.csv / Kidney1MiddleUroSection.csv / Kidney1BottomUroSection.csv
#         Kidney3P_LeftUrothelium.csv / Kidney3P_RightUrothelium.csv
# Output: VisiumHD_kidney_both_deconvolved_collapsed.rds (full object, collapsed
#           RCTD meta.data columns written back — original OBJECT_RDS untouched)
#         DeconvMatrix_Niche9_AllBins.csv / DeconvMatrix_Niche9_MeanBySample.csv
#         Heatmap_Niche9_AllBins_MeanBySample.pdf / Barplot_Niche9_AllBins_MeanBySample.pdf
#         DeconvMatrix_Niche9_FiveRegions.csv / DeconvMatrix_Niche9_MeanByRegion.csv
#         Heatmap_Niche9_FiveRegions_MeanByRegion.pdf / Barplot_Niche9_FiveRegions_MeanByRegion.pdf
################################################################################

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

# ── Paths ─────────────────────────────────────────────────────────────────────
HD_BASE <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSpatialData/VisiumHD"
OUT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/SpatialScripts/output"

OBJECT_RDS          <- file.path(OUT_DIR, "VisiumHD_kidney_both_deconvolved.rds")
OBJECT_RDS_COLLAPSED <- file.path(OUT_DIR, "VisiumHD_kidney_both_deconvolved_collapsed.rds")
NICHE_ID   <- "9"

# ── Collapse related RCTD cell types for interpretability ─────────────────────
VASCULAR_SUBTYPES <- c("rctd_Asc-Vasa-Recta", "rctd_Desc-Vasa-Recta",
                        "rctd_Vas-Efferens", "rctd_Vas-Afferens")

collapse_rctd_celltypes <- function(df) {
  vascular_present <- intersect(VASCULAR_SUBTYPES, colnames(df))
  if (length(vascular_present) > 0) {
    # No na.rm: a bin's rctd_* weights are either all NA (RCTD-excluded bin)
    # or all numeric together, so plain rowSums preserves that NA-ness
    # instead of silently turning excluded bins into a fabricated 0.
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

# ── Write collapsed RCTD columns back onto the object and save ────────────────
# Saved under a new filename rather than overwriting OBJECT_RDS — the RCTD run
# behind that file takes many hours to regenerate, so the uncollapsed 34-celltype
# version is kept intact.
object@meta.data <- meta
message("==> Saving object with collapsed RCTD meta.data ...")
saveRDS(object, OBJECT_RDS_COLLAPSED)
message("  Saved: ", OBJECT_RDS_COLLAPSED)

# ── Helper: heatmap + barplot from a (group x celltype) mean weight matrix ────
plot_composition <- function(mean_mat, group_col_name, out_prefix) {
  mean_df <- as.data.frame(mean_mat)
  mean_df[[group_col_name]] <- rownames(mean_df)
  long_df <- pivot_longer(mean_df, -all_of(group_col_name),
                           names_to = "celltype", values_to = "weight")
  long_df$celltype <- sub("^rctd_", "", long_df$celltype)

  n_groups <- length(unique(long_df[[group_col_name]]))

  # Heatmap: group x cell type, filled by mean RCTD weight
  p_heat <- ggplot(long_df, aes(x = .data[[group_col_name]], y = celltype, fill = weight)) +
    geom_tile(color = "white") +
    scale_fill_gradientn(colors = c("white", "#FED98E", "#FE8929", "#CC4C02", "#7F0000"),
                          name = "Mean\nRCTD weight") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = NULL)
  ggsave(file.path(OUT_DIR, sprintf("Heatmap_%s.pdf", out_prefix)), p_heat,
         width = 4, height = 4)

  # Stacked barplot of cell type composition per group
  p_bar <- ggplot(long_df, aes(x = .data[[group_col_name]], y = weight, fill = celltype)) +
    geom_col(position = "fill") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = "Cell type proportion (mean RCTD weight)", fill = "Cell type")
  ggsave(file.path(OUT_DIR, sprintf("Barplot_%s.pdf", out_prefix)), p_bar,
         width = 5, height = 4)

  message("  Saved: Heatmap_", out_prefix, ".pdf and Barplot_", out_prefix, ".pdf")
}

################################################################################
# Part A: niche 9 across the whole object (both kidneys), grouped by sample_id
################################################################################
message("\n==> Part A: niche 9, whole object, grouped by sample_id ...")

niche9_meta <- meta %>% filter(seurat_cluster.harmony.projected == NICHE_ID)
message(sprintf("  Niche %s bins: %d", NICHE_ID, nrow(niche9_meta)))

niche9_matrix <- niche9_meta[, rctd_cols]
write.csv(niche9_matrix, file.path(OUT_DIR, "DeconvMatrix_Niche9_AllBins.csv"))

mean_by_sample <- niche9_meta %>%
  group_by(sample_id) %>%
  summarise(across(all_of(rctd_cols), ~ mean(.x, na.rm = TRUE))) %>%
  as.data.frame()
rownames(mean_by_sample) <- mean_by_sample$sample_id
mean_by_sample$sample_id <- NULL
write.csv(mean_by_sample, file.path(OUT_DIR, "DeconvMatrix_Niche9_MeanBySample.csv"))

plot_composition(mean_by_sample, "sample_id", "Niche9_AllBins_MeanBySample")

################################################################################
# Part B: niche 9 within the five integrated urothelium regions (OnlyUrotheliumCluster)
################################################################################
message("\n==> Part B: niche 9 within the five integrated urothelium regions ...")

region_specs <- list(
  Pelviscalyceal_urothelium   = list(csv = "Kidney1TopUroSection.csv",     prefix = "kidney_"),
  PelvisCalyxLumen_urothelium = list(csv = "Kidney1MiddleUroSection.csv",  prefix = "kidney_"),
  Periurothelial_urothelium   = list(csv = "Kidney1BottomUroSection.csv",  prefix = "kidney_"),
  ThreeSideLeftUrothelium     = list(csv = "Kidney3P_LeftUrothelium.csv",  prefix = "kidney3p_"),
  ThreeSideRightUrothelium    = list(csv = "Kidney3P_RightUrothelium.csv", prefix = "kidney3p_")
)

region_meta_list <- list()
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
  region_meta <- meta[matched, ]
  region_meta$Region <- region_name
  region_meta_list[[region_name]] <- region_meta
}

five_regions_meta <- do.call(rbind, region_meta_list)

# OnlyUrotheliumCluster: restrict the five integrated regions to niche 9
only_uro_meta <- five_regions_meta %>% filter(seurat_cluster.harmony.projected == NICHE_ID)
message(sprintf("  OnlyUrotheliumCluster (niche %s) bins: %d", NICHE_ID, nrow(only_uro_meta)))
print(table(only_uro_meta$Region))

only_uro_matrix <- only_uro_meta[, rctd_cols]
write.csv(only_uro_matrix, file.path(OUT_DIR, "DeconvMatrix_Niche9_FiveRegions.csv"))

mean_by_region <- only_uro_meta %>%
  group_by(Region) %>%
  summarise(across(all_of(rctd_cols), ~ mean(.x, na.rm = TRUE))) %>%
  as.data.frame()
rownames(mean_by_region) <- mean_by_region$Region
mean_by_region$Region <- NULL
write.csv(mean_by_region, file.path(OUT_DIR, "DeconvMatrix_Niche9_MeanByRegion.csv"))

# Heatmap restricted to cell types with mean weight > 0.01 in at least one region
# (drops rare/noisy RCTD calls that clutter the full heatmap)
celltypes_above_0.01 <- colnames(mean_by_region)[
  apply(mean_by_region, 2, function(x) max(x, na.rm = TRUE) > 0.01)
]
mean_by_region_filtered <- mean_by_region[, celltypes_above_0.01, drop = FALSE]

plot_composition(mean_by_region_filtered, "Region", "Niche9_FiveRegions_MeanByRegion")

message("\n==> Done.")
