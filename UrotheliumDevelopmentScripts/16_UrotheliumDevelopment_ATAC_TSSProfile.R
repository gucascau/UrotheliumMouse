################################################################################
# 16_UrotheliumDevelopment_ATAC_TSSProfile.R
#
# Figure L. Normalized TSS enrichment profile (+/- 1kb around each TSS, i.e.
#   a 2kb window) for each developmental stage, overlaid on one panel.
#
# Figure M. Same profile/panel as Figure L, but with a much wider rolling-
#   mean window (101bp vs. Figure L's 21bp) per the user's explicit ask for
#   "a more flat/smooth line" version -- same underlying data, smoothing
#   only, so each stage reads as a clean peak instead of a jagged line.
#
# Script 08 ran TSSEnrichment(..., fast = TRUE), which only returns the
# scalar per-cell TSS.enrichment score (used already for Figure K) and does
# NOT retain the positional insertion-frequency profile needed to draw an
# actual enrichment curve. Getting the curve requires rerunning
# TSSEnrichment(fast = FALSE), which re-scans Tn5 insertions from the
# fragments file within +/-1000bp of every annotated TSS for all 715 cells --
# confirmed interactively this takes ~10 minutes (not seconds like the rest
# of this pipeline's QC scripts), so this is its own script/job rather than
# folded into 14_UrotheliumDevelopment_ATAC_QC_Figures.R.
#
# TSSEnrichment()/TSSPlot() are both flagged deprecated in this Signac
# version (in favor of ATACqc()) but still fully functional and are what the
# rest of this pipeline (script 08) already standardized on, so they're used
# here too rather than mixing APIs.
#
# positionEnrichment$TSS ends up with 717 rows for 715 cells -- confirmed
# interactively the 2 extra rows are "expected" and "motif" (Signac's own
# internal bookkeeping rows, not real cell barcodes). TSSPlot(group.by="Age")
# is used to build the per-stage profile rather than hand-rolling the
# aggregation, since it already restricts to real cells via the object's
# metadata and matches Signac's own normalization; only its underlying data
# frame ($data: group/position/count/norm.value) is pulled out so the curve
# can be re-themed and re-colored to match this pipeline's STAGE_COLORS
# (Signac's own TSSPlot() facets with its own default palette).
#
# The raw per-bp profile is jagged at this cohort's per-stage cell counts
# (75-159 cells/stage) -- a light centered rolling mean (zoo::rollmean,
# window = 21bp) is applied per stage purely for display; the saved CSV
# keeps the unsmoothed norm.value too.
#
# Input:  UrotheliumDevelopmentScripts/output/UrotheliumOnly_RNA_ATAC.rds
#         (or the cached UrotheliumOnly_RNA_ATAC_TSSprofile.rds below, if it
#         already exists, to skip the ~10-min recompute on reruns)
# Output: UrotheliumOnly_RNA_ATAC_TSSprofile.rds (uro + positionEnrichment,
#           saved so the figures can be restyled later without a 10-min rerun)
#         FigL_UrotheliumDevelopment_ATAC_TSSEnrichmentProfile_ByStage.pdf
#         FigM_UrotheliumDevelopment_ATAC_TSSEnrichmentProfile_Smoothed.pdf
#         UrotheliumDevelopment_ATAC_TSSEnrichmentProfile_ByStage.csv
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(ggplot2)
  library(dplyr)
  library(zoo)
  library(viridisLite)
})

SCRIPT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/UrotheliumDevelopmentScripts"
OUT_DIR    <- file.path(SCRIPT_DIR, "output")

STAGE_ORDER  <- c("E16.5", "P0", "W3", "W12", "W52", "W92")
STAGE_COLORS <- setNames(viridisLite::viridis(length(STAGE_ORDER)), STAGE_ORDER)

# ── Load the combined RNA+ATAC Uro object (cached TSS profile if present) ──
tssprofile_rds <- file.path(OUT_DIR, "UrotheliumOnly_RNA_ATAC_TSSprofile.rds")
if (file.exists(tssprofile_rds)) {
  message("==> Loading cached TSS-profile object (skipping the ~10-min recompute) ...")
  uro <- readRDS(tssprofile_rds)
  uro$Age <- factor(uro$Age, levels = STAGE_ORDER)
} else {
  message("==> Loading Uro RNA+ATAC object ...")
  in_rds <- file.path(OUT_DIR, "UrotheliumOnly_RNA_ATAC.rds")
  if (!file.exists(in_rds)) {
    stop("Missing ", in_rds, " -- run 08_UrotheliumDevelopment_ATAC_Build.R first.")
  }
  uro <- readRDS(in_rds)
  uro$Age <- factor(uro$Age, levels = STAGE_ORDER)
  message(sprintf("  %d cells, %d ATAC peaks", ncol(uro), nrow(uro[["ATAC"]])))

  message("\n==> Running TSSEnrichment(fast = FALSE) (~10 min, scans fragments near every TSS) ...")
  t0 <- Sys.time()
  uro <- TSSEnrichment(uro, assay = "ATAC", fast = FALSE, verbose = TRUE)
  message(sprintf("  Done in %.1f min", as.numeric(Sys.time() - t0, units = "mins")))

  saveRDS(uro, tssprofile_rds)
  message("  Saved: UrotheliumOnly_RNA_ATAC_TSSprofile.rds")
}

################################################################################
# Figure L: normalized TSS enrichment profile by stage
################################################################################
message("\n==> Building Figure L (TSS enrichment profile by stage) ...")

tss_plot_raw <- TSSPlot(uro, assay = "ATAC", group.by = "Age")
profile_df <- tss_plot_raw$data %>%
  rename(Age = group) %>%
  mutate(Age = factor(Age, levels = STAGE_ORDER)) %>%
  arrange(Age, position)

# Light centered smoothing per stage for display only; norm.value (raw) is
# kept alongside in the saved CSV.
profile_df <- profile_df %>%
  group_by(Age) %>%
  mutate(norm.value.smooth = zoo::rollmean(norm.value, k = 21, fill = NA, align = "center")) %>%
  ungroup()

write.csv(profile_df, file.path(OUT_DIR, "UrotheliumDevelopment_ATAC_TSSEnrichmentProfile_ByStage.csv"),
          row.names = FALSE)
message("  Saved: UrotheliumDevelopment_ATAC_TSSEnrichmentProfile_ByStage.csv")

figL <- ggplot(profile_df, aes(x = position, y = norm.value.smooth, color = Age)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = STAGE_COLORS, name = "Stage") +
  labs(x = "Distance from TSS (bp)", y = "Normalized TSS enrichment score",
       title = "Figure L. snATAC-seq TSS enrichment profiles across development",
       subtitle = "+/-1kb window (2kb total) around annotated TSSs; 21bp centered rolling mean") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave(file.path(OUT_DIR, "FigL_UrotheliumDevelopment_ATAC_TSSEnrichmentProfile_ByStage.pdf"),
       figL, width = 7.5, height = 5)
message("  Saved: FigL_UrotheliumDevelopment_ATAC_TSSEnrichmentProfile_ByStage.pdf")

################################################################################
# Figure M: same profile, heavier smoothing (flatter line)
################################################################################
message("\n==> Building Figure M (heavily-smoothed TSS enrichment profile) ...")

SMOOTH_WINDOW_WIDE <- 101
profile_df <- profile_df %>%
  group_by(Age) %>%
  mutate(norm.value.smooth.wide = zoo::rollmean(norm.value, k = SMOOTH_WINDOW_WIDE,
                                                 fill = NA, align = "center")) %>%
  ungroup()

write.csv(profile_df, file.path(OUT_DIR, "UrotheliumDevelopment_ATAC_TSSEnrichmentProfile_ByStage.csv"),
          row.names = FALSE)
message("  Updated: UrotheliumDevelopment_ATAC_TSSEnrichmentProfile_ByStage.csv (added wide-smooth column)")

figM <- ggplot(profile_df, aes(x = position, y = norm.value.smooth.wide, color = Age)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = STAGE_COLORS, name = "Stage") +
  labs(x = "Distance from TSS (bp)", y = "Normalized TSS enrichment score",
       title = "Figure M. snATAC-seq TSS enrichment profiles across development (smoothed)",
       subtitle = sprintf("+/-1kb window (2kb total) around annotated TSSs; %dbp centered rolling mean",
                           SMOOTH_WINDOW_WIDE)) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave(file.path(OUT_DIR, "FigM_UrotheliumDevelopment_ATAC_TSSEnrichmentProfile_Smoothed.pdf"),
       figM, width = 9, height = 5)
message("  Saved: FigM_UrotheliumDevelopment_ATAC_TSSEnrichmentProfile_Smoothed.pdf")

message("\n==> Done.")
