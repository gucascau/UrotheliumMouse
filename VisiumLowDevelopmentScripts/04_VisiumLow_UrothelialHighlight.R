################################################################################
# 04_VisiumLow_UrothelialHighlight.R
#
# Highlight urothelium-containing spots on the Visium tissue images: a
# per-spot urothelial module score (Seurat::AddModuleScore) from the same
# core uroplakin/Krt20 marker panel used as the "Differentiation / barrier"
# set throughout UrotheliumDevelopmentScripts, then a binary highlight
# (score > 0, i.e. above the module score's local background) so urothelial
# spots visually pop out of the tissue regardless of exact score magnitude.
# Visium spots are 55um and capture multiple cells, so this marks spots
# where the urothelial signature is detectable among mixed spot content,
# not necessarily spots that are purely urothelial.
#
# Spatial panels are built manually with ggplot2, same reason as
# 03_VisiumLow_Diagnostic_Plots.R (SpatialDimPlot/SpatialFeaturePlot render a
# blank/gray gradient instead of the tissue image for these objects).
#
# Input:  VisiumLowDevelopmentScripts/output/VisiumLow_harmony_integrated.rds
# Output: VisiumLow_UrothelialHighlight.pdf (score + binary highlight, all 24 samples)
#         VisiumLow_urothelial_spot_counts.csv (per-sample counts/fractions)
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(grid)
  library(viridisLite)
})

SCRIPT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/VisiumLowDevelopmentScripts"
OUT_DIR    <- file.path(SCRIPT_DIR, "output")

STAGE_ORDER <- c("E16.5", "P0", "W3", "W12", "W52", "W92")

# ── Load ─────────────────────────────────────────────────────────────────────
message("==> Loading integrated object ...")
in_rds <- file.path(OUT_DIR, "VisiumLow_harmony_integrated.rds")
if (!file.exists(in_rds)) {
  stop("Missing ", in_rds, " -- run 02_VisiumLow_Integrate_Harmony.R first.")
}
object <- readRDS(in_rds)
object$Age <- factor(object$Age, levels = STAGE_ORDER)

# ── Urothelial module score ─────────────────────────────────────────────────
UROTHELIAL_MARKERS <- c("Upk1a", "Upk1b", "Upk2", "Upk3a", "Upk3b", "Krt20","Krt14","Krt5", "Trp63")
message("\n==> Computing urothelial module score (", paste(UROTHELIAL_MARKERS, collapse = ", "), ") ...")
object <- AddModuleScore(object, features = list(UROTHELIAL_MARKERS),
                          name = "UrothelialScore", assay = "Spatial", seed = 1)
# AddModuleScore appends a trailing "1" (single feature list) -- rename.
colnames(object@meta.data)[colnames(object@meta.data) == "UrothelialScore1"] <- "UrothelialScore"

# Binary highlight threshold: NOT "score > 0". That was tried first and
# rejected -- visual inspection showed the continuous score forms a tight,
# anatomically correct peak hugging the renal pelvis lumen lining, but
# "score > 0" (merely "above the module score's local background control",
# a weak bar) flagged ~15-30% of ALL spots per stage, including large swaths
# of clearly non-urothelial medulla with only marginal signal -- wildly
# inconsistent with the ~0.3-0.8% single-cell urothelial proportion measured
# in the matched snRNA atlas (UrotheliumDevelopmentScripts Fig 3), even
# allowing generously for spot-level dilution (55um Visium spots capture
# multiple mixed cells, so some inflation over the single-cell figure is
# expected, just not 20-40x). Using the top 5% of scores *within each
# sample* instead: adapts to each sample's own score distribution/depth
# rather than an absolute cutoff, and matches the tight peak actually visible
# in the continuous score panel.
UROTHELIAL_TOPQUANTILE <- 0.95
object$UrothelialHighlight <- object@meta.data %>%
  group_by(sample_id) %>%
  mutate(is_uro = UrothelialScore >= quantile(UrothelialScore, UROTHELIAL_TOPQUANTILE)) %>%
  pull(is_uro) %>%
  ifelse("Urothelial", "Other")
object$UrothelialHighlight <- factor(object$UrothelialHighlight, levels = c("Other", "Urothelial"))

# CAVEAT: because the threshold is a fixed top-5% quantile computed
# separately per sample, the resulting pct_urothelial below is ~5% by
# construction for every sample/stage -- it is a visualization aid, NOT a
# comparable measurement of urothelial abundance across stages. Use
# UrotheliumDevelopmentScripts Fig 3 (per-sample Uro proportion from the
# snRNA atlas) for actual cross-stage abundance comparisons.
message("  Urothelial-positive spots per stage (~5% by construction -- not a cross-stage abundance comparison):")
print(object@meta.data %>% group_by(Age) %>%
        summarise(n_total = n(), n_urothelial = sum(UrothelialHighlight == "Urothelial"),
                  pct_urothelial = round(100 * n_urothelial / n_total, 2), .groups = "drop"))

per_sample <- object@meta.data %>%
  group_by(sample_id, Age) %>%
  summarise(n_total = n(), n_urothelial = sum(UrothelialHighlight == "Urothelial"),
            pct_urothelial = round(100 * n_urothelial / n_total, 2), .groups = "drop") %>%
  arrange(Age, sample_id)
write.csv(per_sample, file.path(OUT_DIR, "VisiumLow_urothelial_spot_counts.csv"), row.names = FALSE)
message("\n  Saved: VisiumLow_urothelial_spot_counts.csv")

# ── Manual spatial plot helper (SpatialFeaturePlot/SpatialDimPlot render a
# blank image for these objects -- see 03_VisiumLow_Diagnostic_Plots.R) ─────
plot_spatial_manual <- function(obj, image_name, color_var, title = "",
                                 color_values = NULL, continuous = FALSE,
                                 color_limits = NULL, legend_title = color_var,
                                 pt_size = 1.2) {
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
    theme(plot.background = element_rect(fill = "white", color = NA)) +
    ggtitle(title)

  if (continuous) {
    # color_limits fixed across samples (passed in) so the same color always
    # means the same score value on every page -- without this, each sample
    # auto-scales to its own min/max and "bright yellow" is not comparable
    # from one spatial map to the next.
    p <- p + scale_color_viridis_c(option = "plasma", name = legend_title, limits = color_limits)
  } else if (!is.null(color_values)) {
    p <- p + scale_color_manual(values = color_values, name = legend_title)
  } else {
    p <- p + labs(color = legend_title)
  }
  p
}

HIGHLIGHT_COLORS <- c(Other = "grey85", Urothelial = "#D55E00")  # Okabe-Ito vermillion

# ── UMAP overview: label + highlight urothelial spots among all 68,180 ─────
message("\n==> Building UMAP highlight of urothelial spots (all spots, all samples) ...")
umap_df <- data.frame(Embeddings(object, "umap.harmony"),
                       UrothelialHighlight = object$UrothelialHighlight)
colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
set.seed(1)
umap_df <- umap_df[sample(nrow(umap_df)), ]  # shuffle draw order (see 03's shuffle note)

n_uro_total <- sum(object$UrothelialHighlight == "Urothelial")
centroid <- umap_df[umap_df$UrothelialHighlight == "Urothelial", ] %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))

p_umap_highlight <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = UrothelialHighlight)) +
  geom_point(size = 0.5, alpha = 0.7) +
  scale_color_manual(values = HIGHLIGHT_COLORS, name = NULL) +
  annotate("text", x = centroid$UMAP_1, y = centroid$UMAP_2, label = "Urothelium",
           fontface = "bold", size = 5, color = "black") +
  coord_fixed() +
  theme_classic(base_size = 12) +
  ggtitle(sprintf("Urothelial spots across all samples (top 5%% per sample; %d / %d total, %.1f%%)",
                   n_uro_total, ncol(object), 100 * n_uro_total / ncol(object))) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))

# ── Per-sample highlight plots (score + binary), all 24 samples ────────────
message("\n==> Building spatial highlight plots (24 samples) ...")
# Fixed color scale across all 24 score panels -- global min/max, computed
# once here, so the same plasma color means the same score value on every
# sample's page (otherwise each page auto-scales to its own min/max and
# "bright yellow" isn't comparable from one spatial map to the next).
score_limits <- range(object$UrothelialScore)
message(sprintf("  Urothelial score color scale fixed to [%.3f, %.3f] across all samples",
                score_limits[1], score_limits[2]))

out_pdf <- file.path(OUT_DIR, "VisiumLow_UrothelialHighlight.pdf")
pdf(out_pdf, width = 14, height = 6)

print(p_umap_highlight)

for (samp in unique(object$sample_id)) {
  message("  ", samp)
  sub_obj <- suppressWarnings(subset(object, cells = colnames(object)[object$sample_id == samp]))
  img_name <- Images(sub_obj)[1]
  n_uro <- sum(sub_obj$UrothelialHighlight == "Urothelial")
  pct_uro <- round(100 * n_uro / ncol(sub_obj), 2)

  tryCatch({
    p_score <- plot_spatial_manual(sub_obj, img_name, "UrothelialScore",
                                    title = paste0(samp, " -- urothelial score"),
                                    continuous = TRUE, color_limits = score_limits,
                                    legend_title = "Score", pt_size = 1.2)
    p_highlight <- plot_spatial_manual(sub_obj, img_name, "UrothelialHighlight",
                                        title = sprintf("Urothelial spots: %d / %d (%.2f%%)",
                                                         n_uro, ncol(sub_obj), pct_uro),
                                        color_values = HIGHLIGHT_COLORS,
                                        legend_title = NULL, pt_size = 1.2)
    print(p_score | p_highlight)
  }, error = function(e) {
    message("    Spatial plot failed for ", samp, ": ", e$message)
  })
}

dev.off()
message("\n  Saved: ", out_pdf)
message("\n==> Done.")
