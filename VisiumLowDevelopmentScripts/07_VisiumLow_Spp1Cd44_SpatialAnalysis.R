################################################################################
# 07_VisiumLow_Spp1Cd44_SpatialAnalysis.R
#
# Spatial evidence for the Spp1 (urothelium) / Cd44 (niche) ligand-receptor
# pair across development, on VisiumLow_deconvolved.rds. Unlike the
# UrotheliumDevelopmentScripts/CellCommunicationScripts assay family, this
# object's Spatial assay rownames are already gene symbols (not Ensembl
# IDs) -- confirmed directly (Spp1/Cd44 present, "data" layer already
# log-normalized, distinct from "counts") -- so no symbol-remapping assay
# is needed here, unlike everywhere else in this project.
#
# Focused on 3 user-specified representative sections (not the full 24
# samples, and not a median-fraction auto-pick): E16.5 = GSM8704855,
# P0 = GSM8704850, W12 = GSM8704853 (see FOCUS_SAMPLES below).
#
# Panel 1: Spp1 spatial feature plot at the 3 focus sections -- directly
#   visualizes "Spp1 is acquired postnatally": absent/low at E16.5,
#   induced by P0, sustained at W12.
# Panel 2: Spp1 + Cd44 co-localization at the same 3 sections -- (a)
#   side-by-side Spp1 vs Cd44 spatial plots, and (b) a two-color RGB blend
#   overlay (red=Spp1, green=Cd44, yellow=both) -- shows whether Cd44+
#   spots sit adjacent to/around the Spp1+ urothelial layer, i.e.
#   spatially poised to interact, not just co-expressed in dissociated
#   data.
# Panel 3: Quantification -- Cd44 expression binned by distance from the
#   nearest Uro+ spot (rctd_Uro > URO_THRESHOLD) within each of the 3
#   focus sections. Distance is expressed in "spot units" (each sample's
#   own median nearest-neighbor spot-to-spot distance, K=6, same
#   convention as 06's Moran's I) rather than microns -- the object's
#   VisiumV1 scale factors (spot/hires/lowres/fiducial) don't reliably
#   map to a physical spot_diameter_fullres here (the "spot" and "hires"
#   factors are identical, which is not the standard scalefactors.json
#   layout), so a self-calibrating per-sample unit avoids trusting an
#   ambiguous conversion while still being comparable across
#   samples/stages. Per-stage Spearman correlation (Cd44 vs distance) and
#   a Cd44 ~ dist_units * Age interaction linear model are reported as the
#   statistical backing for the figure.
#
# Requires 05b_VisiumLow_Deconvolution_FixW12W52.R's patch already applied
# (rctd_Uro must have no NAs).
#
# Input:  VisiumLowDevelopmentScripts/output/VisiumLow_deconvolved.rds
# Output: Fig_VisiumLowSpp1Cd44_Panel1_Spp1Expression.pdf
#         Fig_VisiumLowSpp1Cd44_Panel2a_SideBySide.pdf
#         Fig_VisiumLowSpp1Cd44_Panel2b_ColocBlend.pdf
#         Fig_VisiumLowSpp1Cd44_Panel3_DistanceBinnedCd44.pdf
#         VisiumLowSpp1Cd44_per_spot_distance.csv
#         VisiumLowSpp1Cd44_distance_bin_summary.csv
#         VisiumLowSpp1Cd44_distance_correlation_by_stage.csv
#         VisiumLowSpp1Cd44_distance_lm_summary.txt
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(viridisLite)
  library(patchwork)
})

SCRIPT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/VisiumLowDevelopmentScripts"
OUT_DIR    <- file.path(SCRIPT_DIR, "output")

STAGE_ORDER   <- c("E16.5", "P0", "W3", "W12", "W52", "W92")
STAGE_COLORS  <- setNames(viridisLite::viridis(length(STAGE_ORDER)), STAGE_ORDER)
URO_THRESHOLD <- 0.1
KNN_K         <- 6

# User-specified representative sections (not the median-fraction auto-pick).
FOCUS_SAMPLES <- c(
  "E16.5" = "GSM8704855_20210205-NMKE16.5-Fc1U1Z1Bs1",
  "P0"    = "GSM8704850_20210108-NMK0-FM-Fc1U2Z1Bs1",
  "W12"   = "GSM8704853_20210129-AKICT3F-Fc1U1Z1Bs1"
)
FOCUS_STAGES <- names(FOCUS_SAMPLES)

# ── Load ─────────────────────────────────────────────────────────────────────
message("==> Loading VisiumLow deconvolved object ...")
object <- readRDS(file.path(OUT_DIR, "VisiumLow_deconvolved.rds"))
object$Age <- factor(object$Age, levels = STAGE_ORDER)
DefaultAssay(object) <- "Spatial"

if (any(is.na(object$rctd_Uro))) {
  stop("rctd_Uro has NAs -- run 05b_VisiumLow_Deconvolution_FixW12W52.R first.")
}

missing_samples <- setdiff(FOCUS_SAMPLES, unique(object$sample_id))
if (length(missing_samples) > 0) stop("FOCUS_SAMPLES not found in object: ", paste(missing_samples, collapse = ", "))

GENES <- c("Spp1", "Cd44")
missing_genes <- setdiff(GENES, rownames(object[["Spatial"]]))
if (length(missing_genes) > 0) stop("Missing genes: ", paste(missing_genes, collapse = ", "))

expr <- GetAssayData(object, assay = "Spatial", layer = "data")[GENES, ]
object$Spp1_expr <- as.numeric(expr["Spp1", ])
object$Cd44_expr <- as.numeric(expr["Cd44", ])

# ── sample_id -> image map ──────────────────────────────────────────────────
sample_to_image <- setNames(character(length(FOCUS_SAMPLES)), FOCUS_SAMPLES)
for (img_name in Images(object)) {
  tc <- GetTissueCoordinates(object, image = img_name)
  samp <- unique(object$sample_id[rownames(tc)])
  if (samp %in% FOCUS_SAMPLES) sample_to_image[samp] <- img_name
}
rep_lookup <- FOCUS_SAMPLES
message("Focus sections:")
print(data.frame(Age = FOCUS_STAGES, sample_id = FOCUS_SAMPLES, image = sample_to_image[FOCUS_SAMPLES]))

# ── Reusable spatial plotting helpers (pattern from 06) ────────────────────
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
    labs(title = title) +
    theme(plot.title = element_text(hjust = 0.5, size = 11),
          legend.position = "right")
  if (continuous) {
    p <- p + scale_color_viridis_c(option = "plasma", name = legend_title, limits = color_limits)
  } else {
    p <- p + labs(color = legend_title)
  }
  p
}

plot_spatial_identity <- function(obj, image_name, color_hex_var, title = "", pt_size = 1.0) {
  img_s4 <- obj@images[[image_name]]
  coords <- img_s4@coordinates
  coords$color_val <- obj[[color_hex_var, drop = TRUE]][rownames(coords)]
  img_arr <- img_s4@image
  img_dim <- dim(img_arr)
  g <- rasterGrob(img_arr, width = unit(1, "npc"), height = unit(1, "npc"))
  ggplot(coords, aes(x = imagecol, y = imagerow)) +
    annotation_custom(g, xmin = 0, xmax = img_dim[2], ymin = 0, ymax = img_dim[1]) +
    geom_point(aes(color = color_val), size = pt_size) +
    scale_color_identity() +
    scale_y_reverse() +
    coord_fixed(xlim = c(0, img_dim[2]), ylim = c(img_dim[1], 0)) +
    theme_void() +
    labs(title = title) +
    theme(plot.title = element_text(hjust = 0.5, size = 11))
}

# ── Panel 1: Spp1 spatial feature plots across 3 focus sections ─────────────
message("\n==> Panel 1: Spp1 spatial expression across focus sections ...")
rep_cells   <- rownames(object@meta.data)[object$sample_id %in% FOCUS_SAMPLES]
spp1_limits <- range(object$Spp1_expr[rep_cells], na.rm = TRUE)

panel1_plots <- lapply(FOCUS_STAGES, function(st) {
  samp <- rep_lookup[[st]]
  img  <- sample_to_image[[samp]]
  plot_spatial_manual(object, img, "Spp1_expr", title = st,
                       legend_title = "Spp1", color_limits = spp1_limits)
})
panel1 <- wrap_plots(panel1_plots, nrow = 1) +
  plot_annotation(title = "Spp1 expression is acquired postnatally in the urothelium",
                   theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))
ggsave(file.path(OUT_DIR, "Fig_VisiumLowSpp1Cd44_Panel1_Spp1Expression.pdf"), panel1, width = 15, height = 5.5)
message("  Saved Panel 1")

# ── Panel 2: Spp1 + Cd44 co-localization ─────────────────────────────────────
message("\n==> Panel 2: Spp1/Cd44 co-localization ...")
cd44_limits <- range(object$Cd44_expr[rep_cells], na.rm = TRUE)

panel2_side_plots <- unlist(lapply(FOCUS_STAGES, function(st) {
  samp <- rep_lookup[[st]]
  img  <- sample_to_image[[samp]]
  list(
    plot_spatial_manual(object, img, "Spp1_expr", title = paste0(st, " - Spp1"),
                         legend_title = "Spp1", color_limits = spp1_limits),
    plot_spatial_manual(object, img, "Cd44_expr", title = paste0(st, " - Cd44"),
                         legend_title = "Cd44", color_limits = cd44_limits)
  )
}), recursive = FALSE)
panel2_side <- wrap_plots(panel2_side_plots, nrow = 2, byrow = FALSE) +
  plot_annotation(title = "Spp1 (urothelium) and Cd44 (niche) spatial expression",
                   theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))
ggsave(file.path(OUT_DIR, "Fig_VisiumLowSpp1Cd44_Panel2a_SideBySide.pdf"), panel2_side, width = 15, height = 10)

# Two-color RGB blend overlay: R=Spp1, G=Cd44, clipped at 1st/99th pct
# across the 3 focus sections then min-max scaled -- yellow marks overlap.
spp1_q <- quantile(object$Spp1_expr[rep_cells], c(0.01, 0.99), na.rm = TRUE)
cd44_q <- quantile(object$Cd44_expr[rep_cells], c(0.01, 0.99), na.rm = TRUE)
scale01 <- function(x, q) pmin(pmax((x - q[1]) / (q[2] - q[1]), 0), 1)
object$blend_hex <- rgb(
  scale01(object$Spp1_expr, spp1_q),
  scale01(object$Cd44_expr, cd44_q),
  0
)

panel2_blend_plots <- lapply(FOCUS_STAGES, function(st) {
  samp <- rep_lookup[[st]]
  img  <- sample_to_image[[samp]]
  plot_spatial_identity(object, img, "blend_hex", title = st)
})
panel2_blend <- wrap_plots(panel2_blend_plots, nrow = 1) +
  plot_annotation(title = "Spp1 (red) / Cd44 (green) co-localization -- yellow = overlap",
                   theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))
ggsave(file.path(OUT_DIR, "Fig_VisiumLowSpp1Cd44_Panel2b_ColocBlend.pdf"), panel2_blend, width = 15, height = 5.5)
message("  Saved Panel 2 (side-by-side + blend overlay)")

# ── Panel 3: Distance-binned Cd44 expression from the Uro+ boundary ────────
message("\n==> Panel 3: Distance-binned Cd44 quantification ...")

DIST_BREAKS <- c(-Inf, 0, 2, 4, 6, 8, Inf)
DIST_LABELS <- c("Uro+", "0-2", "2-4", "4-6", "6-8", "8+")

compute_distance_df <- function(samp, stage) {
  img <- sample_to_image[[samp]]
  coords <- object@images[[img]]@coordinates[, c("imagecol", "imagerow")]
  cells  <- rownames(coords)
  uro_pos <- object$rctd_Uro[cells] > URO_THRESHOLD
  if (sum(uro_pos) == 0) return(NULL)

  d <- as.matrix(dist(coords))
  diag(d) <- NA

  # Per-sample spot-spacing unit: median K=6 nearest-neighbor distance
  # (same K as 06's Moran's I) -- self-calibrating, avoids relying on the
  # object's ambiguous VisiumV1 scale factors for a physical unit.
  nn_dist <- apply(d, 1, function(x) sort(x, na.last = TRUE)[1:KNN_K])
  unit_dist <- median(nn_dist, na.rm = TRUE)

  dist_to_uro <- apply(d[, uro_pos, drop = FALSE], 1, min, na.rm = TRUE)
  dist_to_uro[uro_pos] <- 0
  dist_units <- dist_to_uro / unit_dist

  data.frame(
    sample_id  = samp,
    Age        = stage,
    cell       = cells,
    dist_units = dist_units,
    is_uro     = uro_pos,
    Cd44_expr  = object$Cd44_expr[cells],
    Spp1_expr  = object$Spp1_expr[cells]
  )
}

dist_df <- bind_rows(lapply(FOCUS_STAGES, function(st) compute_distance_df(rep_lookup[[st]], st)))
dist_df$dist_bin <- cut(ifelse(dist_df$is_uro, -1, dist_df$dist_units),
                         breaks = DIST_BREAKS, labels = DIST_LABELS)
dist_df$Age <- factor(dist_df$Age, levels = FOCUS_STAGES)

write.csv(dist_df, file.path(OUT_DIR, "VisiumLowSpp1Cd44_per_spot_distance.csv"), row.names = FALSE)

bin_summary <- dist_df %>%
  group_by(Age, dist_bin) %>%
  summarise(mean_Cd44 = mean(Cd44_expr), sem_Cd44 = sd(Cd44_expr) / sqrt(n()), n = n(), .groups = "drop")
write.csv(bin_summary, file.path(OUT_DIR, "VisiumLowSpp1Cd44_distance_bin_summary.csv"), row.names = FALSE)

panel3 <- ggplot(bin_summary, aes(x = dist_bin, y = mean_Cd44, color = Age, group = Age)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_Cd44 - sem_Cd44, ymax = mean_Cd44 + sem_Cd44), width = 0.15) +
  scale_color_manual(values = STAGE_COLORS[FOCUS_STAGES]) +
  labs(x = "Distance from Uro+ boundary (spot units)", y = "Mean Cd44 expression (log-norm)",
       title = "Cd44 expression as a function of distance from the urothelial layer",
       color = "Stage") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(file.path(OUT_DIR, "Fig_VisiumLowSpp1Cd44_Panel3_DistanceBinnedCd44.pdf"), panel3, width = 7.5, height = 5.5)
message("  Saved Panel 3")

# ── Statistics ───────────────────────────────────────────────────────────────
message("\n==> Panel 3 statistics ...")
stage_cor <- dist_df %>%
  group_by(Age) %>%
  summarise(spearman_rho = cor(dist_units, Cd44_expr, method = "spearman"),
            p_value = cor.test(dist_units, Cd44_expr, method = "spearman", exact = FALSE)$p.value,
            n = n())
write.csv(stage_cor, file.path(OUT_DIR, "VisiumLowSpp1Cd44_distance_correlation_by_stage.csv"), row.names = FALSE)
message("  Per-stage Spearman correlation (Cd44 vs distance):")
print(stage_cor)

lm_interaction <- lm(Cd44_expr ~ dist_units * Age, data = dist_df)
lm_summary <- summary(lm_interaction)
message("\n  Linear model Cd44 ~ dist_units * Age (interaction test):")
print(lm_summary$coefficients)
capture.output(lm_summary, file = file.path(OUT_DIR, "VisiumLowSpp1Cd44_distance_lm_summary.txt"))

message("\n==> Done.")
