################################################################################
# 21_UrotheliumDevelopment_ATAC_MotifEnrichment.R
#
# Top enriched TF motifs per developmental stage, following Signac's motif
# vignette (https://stuartlab.org/signac/articles/motif_vignette) workflow:
# per-stage differentially-accessible (DA) peaks -> GC-content-matched
# background peak set -> FindMotifs() hypergeometric enrichment test.
#
# This is deliberately NOT a rerun of 10_UrotheliumDevelopment_chromVAR.R.
# chromVAR (10) gives a continuous per-CELL motif-activity z-score for a
# fixed, user-supplied 20-TF panel. This script instead asks an unbiased
# question -- of all 746 JASPAR2020 vertebrate-CORE motifs already scanned
# into uro[["ATAC"]]@motifs by 10's AddMotifs() call, which ones are
# over-represented in each stage's marker peaks -- and does not touch or
# rerun AddMotifs (motifs + RegionStats/GC.percent are already present in
# the chromVAR.rds this loads).
#
# DA peaks: one-vs-rest FindMarkers(test.use="LR", latent.vars="nCount_ATAC")
# per stage, the standard Signac accessibility-DA test (logistic regression
# against a null of read depth, per the vignette) rather than the
# presence/absence proxy used in scripts 14/15/17/20 -- appropriate here
# because FindMotifs' hypergeometric test needs a real foreground/background
# peak split, not a degree-of-specificity ranking.
#
# Background: AccessiblePeaks(uro) (peaks open in >=1 cell anywhere in the
# object, computed once -- constant across stages) subsampled with
# MatchRegionStats() to match the GC-content distribution of each stage's DA
# peaks, exactly as the vignette does. Using a GC-matched background instead
# of the full peak set avoids enrichment artifacts from motifs that are
# simply GC-rich/poor relative to the genome-wide peak GC distribution.
#
# Input:  UrotheliumDevelopmentScripts/output/UrotheliumOnly_RNA_ATAC_chromVAR.rds
# Output: UrotheliumDevelopment_ATAC_MotifEnrichment_AllStages.csv
#         UrotheliumDevelopment_ATAC_MotifEnrichment_TopMotifs.csv
#         FigR_UrotheliumDevelopment_ATAC_MotifEnrichment_Heatmap.pdf
#         FigS_UrotheliumDevelopment_ATAC_MotifEnrichment_TopMotif_Logos.pdf
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(pheatmap)
  library(RColorBrewer)
  library(viridisLite)
})

SCRIPT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/UrotheliumDevelopmentScripts"
OUT_DIR    <- file.path(SCRIPT_DIR, "output")

STAGE_ORDER  <- c("E16.5", "P0", "W3", "W12", "W52", "W92")
STAGE_COLORS <- setNames(viridisLite::viridis(length(STAGE_ORDER)), STAGE_ORDER)

DA_PADJ_MAX          <- 0.05
N_BACKGROUND_MATCH   <- 40000  # vignette default is 50000; capped per-stage to what's available below
TOP_N_MOTIFS_PER_STAGE <- 10

# ── Load combined RNA+ATAC+chromVAR Uro object (motifs already scanned) ────
message("==> Loading combined RNA+ATAC+chromVAR Uro object ...")
in_rds <- file.path(OUT_DIR, "UrotheliumOnly_RNA_ATAC_chromVAR.rds")
if (!file.exists(in_rds)) {
  stop("Missing ", in_rds, " -- run 10_UrotheliumDevelopment_chromVAR.R first.")
}
uro <- readRDS(in_rds)
uro$Age <- factor(uro$Age, levels = STAGE_ORDER)
DefaultAssay(uro) <- "ATAC"
Idents(uro) <- uro$Age

n_motifs <- length(uro[["ATAC"]]@motifs@motif.names)
message(sprintf("  %d peaks x %d cells, %d motifs already scanned (10's AddMotifs)",
                 nrow(uro[["ATAC"]]), ncol(uro[["ATAC"]]), n_motifs))

meta_feature <- uro[["ATAC"]]@meta.features

# ── Background peak universe: peaks open anywhere, constant across stages ─
message("\n==> Building background peak universe (AccessiblePeaks, all stages) ...")
open_peaks <- AccessiblePeaks(uro)
message(sprintf("  %d accessible peaks available as background pool", length(open_peaks)))

################################################################################
# Per-stage: DA peaks (one-vs-rest) -> GC-matched background -> FindMotifs
################################################################################
all_results <- list()

for (st in STAGE_ORDER) {
  message(sprintf("\n==> Stage %s: DA peaks (one-vs-rest, LR test) ...", st))
  da_peaks <- FindMarkers(
    object      = uro,
    ident.1     = st,
    test.use    = "LR",
    latent.vars = "nCount_ATAC",
    only.pos    = TRUE
  )
  top_da_peak <- rownames(da_peaks)[da_peaks$p_val_adj < DA_PADJ_MAX]
  message(sprintf("  %d DA peaks (p_val_adj < %.2g) of %d tested", length(top_da_peak), DA_PADJ_MAX, nrow(da_peaks)))

  if (length(top_da_peak) == 0) {
    message("  No DA peaks survive threshold -- skipping FindMotifs for this stage.")
    next
  }

  peaks_matched <- MatchRegionStats(
    meta.feature  = meta_feature[open_peaks, ],
    query.feature = meta_feature[top_da_peak, ],
    n             = min(N_BACKGROUND_MATCH, length(open_peaks))
  )

  enriched <- FindMotifs(object = uro, features = top_da_peak, background = peaks_matched)
  enriched <- enriched %>%
    mutate(Stage = st, n_DA_peaks = length(top_da_peak)) %>%
    arrange(pvalue)

  message(sprintf("  Top motif: %s (%s), p = %.2e", enriched$motif[1], enriched$motif.name[1], enriched$pvalue[1]))
  all_results[[st]] <- enriched
}

combined <- bind_rows(all_results) %>% mutate(Stage = factor(Stage, levels = STAGE_ORDER))
write.csv(combined, file.path(OUT_DIR, "UrotheliumDevelopment_ATAC_MotifEnrichment_AllStages.csv"), row.names = FALSE)
message(sprintf("\n  Saved: UrotheliumDevelopment_ATAC_MotifEnrichment_AllStages.csv (%d rows)", nrow(combined)))

top_combined <- combined %>%
  group_by(Stage) %>%
  slice_min(pvalue, n = TOP_N_MOTIFS_PER_STAGE, with_ties = FALSE) %>%
  ungroup()
write.csv(top_combined, file.path(OUT_DIR, "UrotheliumDevelopment_ATAC_MotifEnrichment_TopMotifs.csv"), row.names = FALSE)
message(sprintf("  Saved: UrotheliumDevelopment_ATAC_MotifEnrichment_TopMotifs.csv (top %d/stage)", TOP_N_MOTIFS_PER_STAGE))
print(top_combined %>% count(Stage))

################################################################################
# Figure R: -log10(p) heatmap, union of each stage's top motifs across all
# six stages (every motif was tested in every stage's FindMotifs() call, so
# no NAs -- rows ordered by which stage each motif peaks in, as in Figure C)
################################################################################
if (nrow(top_combined) == 0) {
  message("\n  No stage produced any enriched motif -- skipping Figure R/S.")
} else {
  message("\n==> Building Figure R (motif enrichment heatmap) ...")

  union_motifs <- unique(top_combined$motif)
  motif_name_lookup <- setNames(combined$motif.name, combined$motif)

  neglog10_mat <- sapply(STAGE_ORDER, function(st) {
    df <- all_results[[st]]
    if (is.null(df)) return(rep(NA_real_, length(union_motifs)))
    -log10(df$pvalue[match(union_motifs, df$motif)])
  })
  rownames(neglog10_mat) <- make.unique(motif_name_lookup[union_motifs])
  colnames(neglog10_mat) <- STAGE_ORDER

  finite_max <- max(neglog10_mat[is.finite(neglog10_mat)])
  neglog10_mat[!is.finite(neglog10_mat)] <- finite_max  # cap pvalue == 0

  row_order <- order(apply(neglog10_mat, 1, which.max))
  neglog10_mat <- neglog10_mat[row_order, , drop = FALSE]

  row_anno <- data.frame(top_in_stage = STAGE_ORDER[apply(neglog10_mat, 1, which.max)],
                          row.names = rownames(neglog10_mat))
  anno_colors <- list(top_in_stage = STAGE_COLORS)

  hm_palette <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)

  pdf(file.path(OUT_DIR, "FigR_UrotheliumDevelopment_ATAC_MotifEnrichment_Heatmap.pdf"),
      width = 7, height = max(4, 0.25 * nrow(neglog10_mat)))
  pheatmap(neglog10_mat, cluster_rows = FALSE, cluster_cols = FALSE,
           color = hm_palette, border_color = NA,
           annotation_row = row_anno, annotation_colors = anno_colors,
           fontsize_row = 8,
           main = sprintf("Figure R. Top motif enrichment by stage (-log10 p, top %d/stage, n=%d motifs)",
                           TOP_N_MOTIFS_PER_STAGE, nrow(neglog10_mat)))
  dev.off()
  message("  Saved: FigR_UrotheliumDevelopment_ATAC_MotifEnrichment_Heatmap.pdf")

  ##############################################################################
  # Figure S: sequence logos for the single top motif per stage.
  #
  # Built directly from GetMotifData(slot="pwm") + ggseqlogo rather than via
  # Signac::MotifPlot(): MotifPlot unconditionally overwrites whatever names
  # its `motifs` argument carries with GetMotifData(slot="motif.names")
  # (confirmed in Signac 1.17.1 source), discarding any stage label we'd
  # attach -- and if the same TF motif happens to top two different stages,
  # the two logo panels would collapse under one identical facet name.
  # Prefixing each entry's list name with its Stage (done manually here,
  # ggseqlogo facets by list name) keeps every panel distinct and labeled.
  ##############################################################################
  message("\n==> Building Figure S (top motif sequence logos per stage) ...")
  top1_per_stage <- combined %>%
    group_by(Stage) %>%
    slice_min(pvalue, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(Stage)

  pwm_list <- GetMotifData(uro, assay = "ATAC", slot = "pwm")[top1_per_stage$motif]
  names(pwm_list) <- sprintf("%s: %s", top1_per_stage$Stage, top1_per_stage$motif.name)

  pdf(file.path(OUT_DIR, "FigS_UrotheliumDevelopment_ATAC_MotifEnrichment_TopMotif_Logos.pdf"),
      width = 10, height = 6)
  print(ggseqlogo::ggseqlogo(pwm_list) +
          ggplot2::labs(title = "Figure S. Top enriched motif per developmental stage"))
  dev.off()
  message("  Saved: FigS_UrotheliumDevelopment_ATAC_MotifEnrichment_TopMotif_Logos.pdf")
}

message("\n==> Done.")
