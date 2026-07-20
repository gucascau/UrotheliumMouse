################################################################################
# 06_VisiumLow_UrothelialSpatialAnalyses.R
#
# Four spatial analyses of urothelial RCTD deconvolution across the
# developmental Visium atlas (VisiumLow_deconvolved.rds, all 6 stages --
# requires 05b_VisiumLow_Deconvolution_FixW12W52.R's patch to have been
# applied first; that script fixed a bug that left W12/W52 with NA rctd_*
# values in the original 05 run).
#
# Panel A: Representative H&E section per stage, colored by rctd_Uro
#   proportion (continuous), not a module score. "Representative" = the
#   sample whose per-sample Uro+ spot fraction (Panel B's core statistic) is
#   closest to that stage's median across replicates -- an objective pick,
#   not an arbitrary one; the chosen sample_id is in each panel's title.
# Panel B: Per-kidney Uro+ spot fraction across stage (bar = mean +/- SEM,
#   points = individual kidneys), same style as UrotheliumDevelopmentScripts/
#   01_UrotheliumDevelopment_Figures.R's Figure 3. Kruskal-Wallis test for a
#   stage effect (nonparametric -- unbalanced small groups, n=2-10/stage).
#   "Uro+" spot = rctd_Uro > URO_THRESHOLD (0.1); mean Uro proportion among
#   Uro+ spots plotted alongside as a secondary panel.
# Panel C: Per-kidney spatial consolidation index across stage -- Moran's I
#   of rctd_Uro proportion over a k=6 nearest-neighbor spatial graph (k=6 --
#   a Visium hex grid's neighbor count). High positive Moran's I = spots
#   with similar Uro proportion cluster together (consolidated, contiguous
#   structure); near-zero/negative = scattered/diffuse. Meant to sit next to
#   the pseudotime plateau figure (UrotheliumDevelopmentScripts Fig 2) on
#   the same x-axis (stage) -- not combined into one figure here, just
#   designed to be read side by side.
# Panel D: Mean cell-type composition of Uro+ spots per stage (stacked bar,
#   includes Uro itself + top co-occurring types + "Other"). Spatial
#   evidence (independent of the dissociated scRNA-seq marker analysis) for
#   whether urothelium-containing spots co-mix with nephric-progenitor-like
#   signal (NP/UBP/IM) early and shift to mature/stromal signal later.
#   CAVEAT: per 05_VisiumLow_Deconvolution.R's MIN_CELLS_PER_TYPE filter,
#   rare cell types (e.g. UBP, IM_proliferate, LOH_AL_proliferating) are
#   dropped from whichever stage-specific reference didn't have >= 25 cells
#   of that type -- a "0" for such a type at some stage can mean "excluded,
#   too rare to model" rather than "genuinely absent." Not silently hidden:
#   flagged in this script's printed output before Panel D is built.
#
# Input:  VisiumLowDevelopmentScripts/output/VisiumLow_deconvolved.rds
# Output: Fig_VisiumLowUro_A_RepresentativeSections.pdf
#         Fig_VisiumLowUro_B_SpotFractionByStage.pdf
#         Fig_VisiumLowUro_C_SpatialConsolidationIndex.pdf
#         Fig_VisiumLowUro_D_CoLocalizationComposition.pdf
#         VisiumLowUro_per_sample_stats.csv (per-kidney stats used in B/C)
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(grid)
  library(ape)
  library(viridisLite)
})

SCRIPT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/VisiumLowDevelopmentScripts"
OUT_DIR    <- file.path(SCRIPT_DIR, "output")

STAGE_ORDER   <- c("E16.5", "P0", "W3", "W12", "W52", "W92")
STAGE_COLORS  <- setNames(viridisLite::viridis(length(STAGE_ORDER)), STAGE_ORDER)
URO_THRESHOLD <- 0.1
KNN_K         <- 6

message("==> Loading deconvolved object ...")
in_rds <- file.path(OUT_DIR, "VisiumLow_deconvolved.rds")
if (!file.exists(in_rds)) stop("Missing ", in_rds)
object <- readRDS(in_rds)
object$Age <- factor(object$Age, levels = STAGE_ORDER)
if (any(is.na(object$rctd_Uro))) {
  stop("rctd_Uro still has NAs -- run 05b_VisiumLow_Deconvolution_FixW12W52.R first.")
}
message(sprintf("  %d spots, %d samples", ncol(object), length(unique(object$sample_id))))

sample_age_map <- object@meta.data %>% distinct(sample_id, Age)
stopifnot(!any(duplicated(sample_age_map$sample_id)))

# Map sample_id -> image name without ever calling subset() on the full
# (large) merged object -- each image's @coordinates already only contains
# that sample's own spots, so plot_spatial_manual() below can index
# `object@images[[img_name]]` directly on the full object.
sample_to_image <- setNames(character(nrow(sample_age_map)), sample_age_map$sample_id)
for (img_name in Images(object)) {
  tc <- GetTissueCoordinates(object, image = img_name)
  samp <- unique(object$sample_id[rownames(tc)])
  stopifnot(length(samp) == 1)
  sample_to_image[samp] <- img_name
}

# ── Manual spatial plot helper (SpatialDimPlot renders blank for these ──────
# ── objects -- same issue/fix documented in 03_VisiumLow_Diagnostic_Plots.R) ─
plot_spatial_manual <- function(obj, image_name, color_var, title = "",
                                 continuous = TRUE, legend_title = color_var,
                                 pt_size = 1.0, color_limits = NULL) {
  img_s4 <- obj@images[[image_name]]
  coords <- img_s4@coordinates
  coords$color_val <- obj[[color_var, drop = TRUE]][rownames(coords)]
  img_arr <- img_s4@image
  img_dim <- dim(img_arr)
  g <- rasterGrob(img_arr, width = unit(1, "npc"), height = unit(1, "npc"))

  p <- ggplot(coords, aes(x = imagecol, y = imagerow)) +
    annotation_custom(g, xmin = 0, xmax = img_dim[2], ymin = 0, ymax = img_dim[1]) +
    geom_point(aes(color = color_val), size = pt_size) +
    scale_y_reverse() +
    coord_fixed(xlim = c(0, img_dim[2]), ylim = c(img_dim[1], 0)) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA),
          plot.title = element_text(size = 9, hjust = 0.5)) +
    ggtitle(title)

  if (continuous) {
    p <- p + scale_color_viridis_c(option = "plasma", name = legend_title, limits = color_limits)
  } else {
    p <- p + labs(color = legend_title)
  }
  p
}

# ── Per-sample statistics used across Panels A-C ────────────────────────────
message("\n==> Computing per-sample statistics ...")
per_sample_stats <- object@meta.data %>%
  group_by(sample_id, Age) %>%
  summarise(
    n_spots            = n(),
    n_uro_pos          = sum(rctd_Uro > URO_THRESHOLD),
    frac_uro_pos       = n_uro_pos / n_spots,
    mean_uro_among_pos = mean(rctd_Uro[rctd_Uro > URO_THRESHOLD]),
    .groups = "drop"
  ) %>%
  mutate(mean_uro_among_pos = ifelse(is.nan(mean_uro_among_pos), NA_real_, mean_uro_among_pos))

# ── Moran's I spatial consolidation index per sample (feeds Panel C) ───────
message("==> Computing Moran's I (spatial consolidation) per sample ...")
compute_morans_i <- function(samp) {
  img_name <- sample_to_image[[samp]]
  coords <- object@images[[img_name]]@coordinates[, c("imagecol", "imagerow")]
  vals <- object$rctd_Uro[rownames(coords)]

  d <- as.matrix(dist(coords))
  n <- nrow(coords)
  w <- matrix(0, nrow = n, ncol = n)
  for (i in seq_len(n)) {
    nn <- order(d[i, ])[2:(KNN_K + 1)]  # rank 1 is self (distance 0)
    w[i, nn] <- 1
  }
  w <- pmax(w, t(w))  # symmetrize: a neighbor relationship counts either direction

  mi <- tryCatch(ape::Moran.I(vals, w),
                 error = function(e) list(observed = NA_real_, p.value = NA_real_))
  data.frame(sample_id = samp, morans_i = mi$observed, morans_i_p = mi$p.value)
}
morans_df <- bind_rows(lapply(names(sample_to_image), compute_morans_i))
per_sample_stats <- per_sample_stats %>% left_join(morans_df, by = "sample_id")

out_csv <- file.path(OUT_DIR, "VisiumLowUro_per_sample_stats.csv")
write.csv(per_sample_stats, out_csv, row.names = FALSE)
message("  Saved: ", out_csv)
message("\n  Per-sample statistics:")
print(as.data.frame(per_sample_stats %>% arrange(Age)))

################################################################################
# Panel A: Representative section per stage
################################################################################
message("\n==> Panel A: representative section per stage ...")

representative_samples <- per_sample_stats %>%
  group_by(Age) %>%
  mutate(stage_median = median(frac_uro_pos)) %>%
  slice_min(abs(frac_uro_pos - stage_median), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(Age)
message("  Representative sample chosen per stage (closest to stage median Uro+ fraction):")
print(as.data.frame(representative_samples %>% select(Age, sample_id, frac_uro_pos, stage_median)))

uro_range <- range(object$rctd_Uro)  # shared color scale across all 6 panels

panelA_plots <- lapply(seq_len(nrow(representative_samples)), function(i) {
  samp <- representative_samples$sample_id[i]
  img_name <- sample_to_image[[samp]]
  plot_spatial_manual(
    object, img_name, "rctd_Uro",
    title = sprintf("%s  (%s)\nUro+ frac = %.1f%%",
                     representative_samples$Age[i], samp,
                     100 * representative_samples$frac_uro_pos[i]),
    continuous = TRUE, legend_title = "RCTD\nUro\nproportion",
    pt_size = 1.0, color_limits = uro_range
  )
})
figA <- wrap_plots(panelA_plots, nrow = 2) +
  plot_annotation(
    title = "Panel A. Representative section per stage, colored by deconvolved urothelial proportion",
    subtitle = "Representative = sample closest to that stage's median Uro+ spot fraction (sample_id noted per panel)"
  )
ggsave(file.path(OUT_DIR, "Fig_VisiumLowUro_A_RepresentativeSections.pdf"), figA, width = 16, height = 10)
message("  Saved: Fig_VisiumLowUro_A_RepresentativeSections.pdf")

# ── Supplementary: all 24 samples, same fixed color scale ───────────────────
# Every kidney (not just the one representative pick per stage), all on the
# same rctd_Uro color scale (uro_range, shared with the 6-panel figure above)
# so replicates are directly comparable to each other and to the
# representative pick -- a QC/completeness companion to the curated 6-panel
# figure, not a replacement for it.
message("\n==> Panel A (supplementary): all 24 samples, shared color scale ...")
all_samples_ordered <- per_sample_stats %>% arrange(Age, sample_id)

panelA_all_plots <- lapply(seq_len(nrow(all_samples_ordered)), function(i) {
  samp <- all_samples_ordered$sample_id[i]
  img_name <- sample_to_image[[samp]]
  plot_spatial_manual(
    object, img_name, "rctd_Uro",
    title = sprintf("%s\n%s", all_samples_ordered$Age[i], samp),
    continuous = TRUE, legend_title = "RCTD\nUro\nproportion",
    pt_size = 0.6, color_limits = uro_range
  )
})
figA_all <- wrap_plots(panelA_all_plots, ncol = 5) +
  plot_annotation(
    title = "Panel A (all samples). All 24 kidneys, colored by deconvolved urothelial proportion",
    subtitle = "Same fixed RCTD Uro-proportion color scale across every sample for direct comparison"
  )
ggsave(file.path(OUT_DIR, "Fig_VisiumLowUro_A_AllSamples.pdf"), figA_all, width = 20, height = 24, limitsize = FALSE)
message("  Saved: Fig_VisiumLowUro_A_AllSamples.pdf")

################################################################################
# Panel B: Quantified urothelial spot fraction across stages
################################################################################
message("\n==> Panel B: quantified urothelial spot fraction across stages ...")

kw_frac <- kruskal.test(frac_uro_pos ~ Age, data = per_sample_stats)
message(sprintf("  Kruskal-Wallis (Uro+ spot fraction ~ stage): chi-sq=%.2f, df=%d, p=%.4g",
                 kw_frac$statistic, kw_frac$parameter, kw_frac$p.value))

stage_summary_B1 <- per_sample_stats %>%
  group_by(Age) %>%
  summarise(mean_frac = mean(frac_uro_pos), sem_frac = sd(frac_uro_pos) / sqrt(n()), .groups = "drop")

panelB1 <- ggplot(stage_summary_B1, aes(x = Age, y = 100 * mean_frac, fill = Age)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = 100 * (mean_frac - sem_frac), ymax = 100 * (mean_frac + sem_frac)),
                width = 0.15, linewidth = 0.4) +
  geom_jitter(data = per_sample_stats, aes(x = Age, y = 100 * frac_uro_pos),
              inherit.aes = FALSE, width = 0.1, size = 1.8, alpha = 0.7) +
  scale_fill_manual(values = STAGE_COLORS, guide = "none") +
  labs(x = "Developmental stage", y = sprintf("Uro+ spots (%% of spots, rctd_Uro > %.2f)", URO_THRESHOLD),
       title = "Panel B1. Urothelial spot fraction across development",
       subtitle = sprintf("Bars = mean +/- SEM across kidneys; points = individual kidneys. Kruskal-Wallis p = %.3g", kw_frac$p.value)) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13))

per_sample_valid <- per_sample_stats %>% filter(!is.na(mean_uro_among_pos))
kw_mean <- kruskal.test(mean_uro_among_pos ~ Age, data = per_sample_valid)
stage_summary_B2 <- per_sample_valid %>%
  group_by(Age) %>%
  summarise(mean_val = mean(mean_uro_among_pos), sem_val = sd(mean_uro_among_pos) / sqrt(n()), .groups = "drop")

panelB2 <- ggplot(stage_summary_B2, aes(x = Age, y = mean_val, fill = Age)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean_val - sem_val, ymax = mean_val + sem_val), width = 0.15, linewidth = 0.4) +
  geom_jitter(data = per_sample_valid, aes(x = Age, y = mean_uro_among_pos),
              inherit.aes = FALSE, width = 0.1, size = 1.8, alpha = 0.7) +
  scale_fill_manual(values = STAGE_COLORS, guide = "none") +
  labs(x = "Developmental stage", y = "Mean RCTD Uro proportion (among Uro+ spots)",
       title = "Panel B2. Urothelial signal strength among Uro+ spots",
       subtitle = sprintf("Kruskal-Wallis p = %.3g", kw_mean$p.value)) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13))

figB <- panelB1 / panelB2
ggsave(file.path(OUT_DIR, "Fig_VisiumLowUro_B_SpotFractionByStage.pdf"), figB, width = 7, height = 9)
message("  Saved: Fig_VisiumLowUro_B_SpotFractionByStage.pdf")

################################################################################
# Panel C: Spatial dispersion/consolidation index across stages
################################################################################
message("\n==> Panel C: spatial consolidation index (Moran's I) across stages ...")

kw_mi <- kruskal.test(morans_i ~ Age, data = per_sample_stats)
message(sprintf("  Kruskal-Wallis (Moran's I ~ stage): chi-sq=%.2f, df=%d, p=%.4g",
                 kw_mi$statistic, kw_mi$parameter, kw_mi$p.value))

figC <- ggplot(per_sample_stats, aes(x = Age, y = morans_i)) +
  geom_violin(aes(fill = Age), alpha = 0.5, color = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +
  scale_fill_manual(values = STAGE_COLORS, guide = "none") +
  labs(x = "Developmental stage", y = "Moran's I (spatial autocorrelation of rctd_Uro)",
       title = "Panel C. Spatial consolidation of urothelial signal across development",
       subtitle = sprintf(
         "Higher Moran's I = urothelial spots spatially clustered (consolidated); ~0 = scattered (diffuse). k=%d nearest-neighbor graph. Kruskal-Wallis p = %.3g",
         KNN_K, kw_mi$p.value)) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13))
ggsave(file.path(OUT_DIR, "Fig_VisiumLowUro_C_SpatialConsolidationIndex.pdf"), figC, width = 7, height = 5)
message("  Saved: Fig_VisiumLowUro_C_SpatialConsolidationIndex.pdf")

################################################################################
# Panel D: Co-localization / neighboring composition
################################################################################
message("\n==> Panel D: co-localization composition of Uro+ spots ...")

rctd_cols <- setdiff(grep("^rctd_", colnames(object@meta.data), value = TRUE), "rctd_dominant_celltype")

uro_pos_meta <- object@meta.data %>%
  filter(rctd_Uro > URO_THRESHOLD) %>%
  select(Age, all_of(rctd_cols))

# Cell types dropped from at least one stage's reference (structural 0 there,
# not necessarily biologically absent) -- flagged before Panel D, not hidden.
zero_by_stage <- uro_pos_meta %>%
  group_by(Age) %>%
  summarise(across(all_of(rctd_cols), ~ all(. == 0)), .groups = "drop")
message("  Cell types with structural all-zero at some (not all) stages among Uro+ spots")
message("  (per 05's MIN_CELLS_PER_TYPE reference filter -- '0' there means 'excluded from")
message("  that stage's reference', not necessarily 'absent from the tissue'):")
for (ct in rctd_cols) {
  flagged_stages <- as.character(zero_by_stage$Age[zero_by_stage[[ct]]])
  if (length(flagged_stages) > 0 && length(flagged_stages) < length(STAGE_ORDER)) {
    message(sprintf("    %s: all-zero at %s", sub("^rctd_", "", ct), paste(flagged_stages, collapse = ", ")))
  }
}

comp_by_stage <- uro_pos_meta %>%
  group_by(Age) %>%
  summarise(across(all_of(rctd_cols), mean), .groups = "drop") %>%
  pivot_longer(-Age, names_to = "celltype", values_to = "mean_proportion") %>%
  mutate(celltype = sub("^rctd_", "", celltype))

# Keep the top types by overall mean proportion, always keeping the
# nephric-progenitor axis (Uro itself + NP/UBP/IM family) the user is
# specifically asking about even if individually small; everything else -> "Other".
overall_rank <- comp_by_stage %>% group_by(celltype) %>% summarise(m = mean(mean_proportion), .groups = "drop") %>% arrange(desc(m))
always_keep <- intersect(c("Uro", "NP", "NP_proliferate", "UBP", "IM", "IM_proliferate"), unique(comp_by_stage$celltype))
top_types <- union(always_keep, head(overall_rank$celltype, 8))

comp_by_stage <- comp_by_stage %>%
  mutate(celltype_grouped = ifelse(celltype %in% top_types, celltype, "Other")) %>%
  group_by(Age, celltype_grouped) %>%
  summarise(mean_proportion = sum(mean_proportion), .groups = "drop")

plot_types <- unique(comp_by_stage$celltype_grouped)
non_other <- setdiff(plot_types, "Other")
fill_colors <- setNames(c(viridisLite::viridis(length(non_other), option = "turbo"), "grey80"),
                         c(non_other, "Other"))

figD <- ggplot(comp_by_stage, aes(x = Age, y = mean_proportion, fill = celltype_grouped)) +
  geom_col(position = "stack", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = fill_colors, name = "Cell type") +
  labs(x = "Developmental stage", y = "Mean RCTD proportion (Uro+ spots)",
       title = "Panel D. Co-localization composition of urothelial (Uro+) spots",
       subtitle = sprintf(
         "Uro+ = rctd_Uro > %.2f. Includes Uro itself; nephric-progenitor types (NP/NP_proliferate/UBP/IM/IM_proliferate) always broken out when present.",
         URO_THRESHOLD)) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13))
ggsave(file.path(OUT_DIR, "Fig_VisiumLowUro_D_CoLocalizationComposition.pdf"), figD, width = 8, height = 6)
message("  Saved: Fig_VisiumLowUro_D_CoLocalizationComposition.pdf")

message("\n==> Done.")
