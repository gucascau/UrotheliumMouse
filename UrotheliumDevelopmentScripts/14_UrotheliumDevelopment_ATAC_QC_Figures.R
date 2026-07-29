################################################################################
# 14_UrotheliumDevelopment_ATAC_QC_Figures.R
#
# Three snATAC-seq QC/summary figures for the combined RNA+ATAC Uro object
# built by 08_UrotheliumDevelopment_ATAC_Build.R.
#
# Figure H. Per-cell fragment-count distribution (log10 nCount_ATAC) as a
#   violin by developmental stage -- the standard ATAC library-complexity QC
#   view. nCount_ATAC/atac_total_counts is colSums() of the peak x cell
#   counts matrix (fragments-in-peaks per cell, computed once in script 08);
#   the two columns are identical, nCount_ATAC is used since it's Signac's
#   own name for this metric.
#
# Figure I. Peak-set genomic annotation, stacked by developmental stage
#   (Promoter / Exonic / Intronic / Distal segments per stage bar), using
#   ChIPseeker's standard priority order (Promoter > Exon > Intron >
#   Downstream/Intergenic).
#
#   ChIPseeker::annotatePeak() is called with TxDb = EnsDb.Mmusculus.v79
#   directly -- ensembldb registers the GenomicFeatures generics annotatePeak
#   needs (transcripts/exonsBy/etc.), so no separate TxDb.Mmusculus.UCSC.*
#   package is required (confirmed interactively; that package isn't
#   installed under module R/4.4.0 here, but ChIPseeker + EnsDb.Mmusculus.v79
#   are). The peak GRanges carry UCSC-style seqnames ("chr19") from script 08
#   while EnsDb.Mmusculus.v79 is Ensembl-style ("19") -- annotatePeak() does
#   NOT reconcile this itself (silently returns 0 overlaps / "invalid
#   subscript" errors); seqlevelsStyle() is forced to "Ensembl" on a copy of
#   the peaks before annotating, checked against a 3-peak smoke test first.
#
#   Per-stage peak membership: the original MACS peak_called_in field (set by
#   CallPeaks(group.by="Age") in script 08) does NOT survive into this RDS --
#   FeatureMatrix()/CreateChromatinAssay() rebuild the ChromatinAssay from a
#   plain counts matrix with no mcols, so granges(uro[["ATAC"]]) here carries
#   zero metadata columns (confirmed interactively: peak_called_in grepl
#   matched 0/29310 peaks for every stage), and regenerating the true MACS
#   call would mean rerunning script 08's ~45-60 min CallPeaks job. Per the
#   user's choice, "peak present in stage X" is instead read off the counts
#   matrix already in this object: a peak counts toward a stage if it has
#   nonzero accessibility in >= MIN_CELLS_DETECTED cells of that stage
#   (Signac's own CreateChromatinAssay default min.cells threshold). This is
#   a coverage-based proxy for the pipeline's original per-stage MACS calls,
#   not identical to them, but needs no rerun. A peak can be "present" in
#   more than one stage (same non-exclusive convention CallPeaks used), which
#   is fine here since each stage gets its own independent stacked bar.
#
# Figure K. Per-cell TSS enrichment score (uro$TSS.enrichment) distribution by
#   developmental stage, same violin+box style as Figure H. TSS.enrichment
#   was already computed once in script 08 via TSSEnrichment(fast=TRUE) (the
#   scalar per-cell score only, no positional profile) -- reused as-is here,
#   not recomputed.
#
# Input:  UrotheliumDevelopmentScripts/output/UrotheliumOnly_RNA_ATAC.rds
# Output: FigH_UrotheliumDevelopment_ATAC_FragmentCounts_ViolinByStage.pdf
#         FigI_UrotheliumDevelopment_ATAC_PeakAnnotation_Barplot.pdf
#         FigK_UrotheliumDevelopment_ATAC_TSSEnrichmentScore_ViolinByStage.pdf
#         UrotheliumDevelopment_ATAC_PeakAnnotation_Summary.csv
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(EnsDb.Mmusculus.v79)
  library(ChIPseeker)
  library(ggplot2)
  library(dplyr)
  library(viridisLite)
})

SCRIPT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/UrotheliumDevelopmentScripts"
OUT_DIR    <- file.path(SCRIPT_DIR, "output")

STAGE_ORDER  <- c("E16.5", "P0", "W3", "W12", "W52", "W92")
STAGE_COLORS <- setNames(viridisLite::viridis(length(STAGE_ORDER)), STAGE_ORDER)

# ── Load the combined RNA+ATAC Uro object ───────────────────────────────────
message("==> Loading Uro RNA+ATAC object ...")
in_rds <- file.path(OUT_DIR, "UrotheliumOnly_RNA_ATAC.rds")
if (!file.exists(in_rds)) {
  stop("Missing ", in_rds, " -- run 08_UrotheliumDevelopment_ATAC_Build.R first.")
}
uro <- readRDS(in_rds)
uro$Age <- factor(uro$Age, levels = STAGE_ORDER)
message(sprintf("  %d cells, %d ATAC peaks", ncol(uro), nrow(uro[["ATAC"]])))

################################################################################
# Figure H: snATAC-seq fragment-count distribution by developmental stage
################################################################################
message("\n==> Building Figure H (fragment-count violin by stage) ...")

qc_df <- data.frame(Age = uro$Age, nCount_ATAC = uro$nCount_ATAC) %>%
  mutate(log10_nCount_ATAC = log10(nCount_ATAC))

message("  nCount_ATAC summary by stage:")
print(qc_df %>% group_by(Age) %>%
        summarise(median = median(nCount_ATAC), n = n(), .groups = "drop"))

figH <- ggplot(qc_df, aes(x = Age, y = log10_nCount_ATAC, fill = Age)) +
  geom_violin(trim = TRUE, linewidth = 0.3, color = "black") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, linewidth = 0.3) +
  scale_fill_manual(values = STAGE_COLORS, guide = "none") +
  labs(x = "Developmental stage", y = expression(log[10]~"(ATAC fragment counts)"),
       title = "Figure H. snATAC-seq fragment counts across development",
       subtitle = "Per-cell fragments-in-peaks (nCount_ATAC); box = median/IQR") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave(file.path(OUT_DIR, "FigH_UrotheliumDevelopment_ATAC_FragmentCounts_ViolinByStage.pdf"),
       figH, width = 6.5, height = 4.5)
message("  Saved: FigH_UrotheliumDevelopment_ATAC_FragmentCounts_ViolinByStage.pdf")

################################################################################
# Figure I: Peak-set genomic annotation, stacked by developmental stage
################################################################################
message("\n==> Building Figure I (peak genomic annotation, stacked by stage) ...")

peaks <- granges(uro[["ATAC"]])
peaks_ens <- peaks
seqlevelsStyle(peaks_ens) <- "Ensembl"
message(sprintf("  Annotating %d merged peaks ...", length(peaks_ens)))

anno <- annotatePeak(peaks_ens, TxDb = EnsDb.Mmusculus.v79,
                      tssRegion = c(-3000, 3000), verbose = FALSE)
anno_df <- as.data.frame(anno)

# Collapse ChIPseeker's detailed Feature levels into the 4 requested buckets.
# UTRs are exonic sequence -> Exonic; Downstream (near 3' gene end, ChIPseeker's
# own catch-all alongside Distal Intergenic) -> Distal, same as the paper
# convention of calling everything outside promoter/gene body "distal".
collapse_feature <- function(x) {
  dplyr::case_when(
    grepl("^Promoter", x)        ~ "Promoter",
    grepl("Exon|UTR", x)         ~ "Exonic",
    grepl("Intron", x)           ~ "Intronic",
    TRUE                          ~ "Distal"
  )
}
ANNO_LEVELS <- c("Promoter", "Exonic", "Intronic", "Distal")
peak_category <- collapse_feature(anno_df$annotation)  # length == length(peaks), same order

# Okabe-Ito colorblind-safe categorical set, fixed order (consistent with
# CELLCLASS_COLORS in 01_UrotheliumDevelopment_Figures.R).
ANNO_COLORS <- setNames(c("#D55E00", "#009E73", "#0072B2", "#B3B3B3"), ANNO_LEVELS)

# ── Which peaks are "present" in each stage (counts-matrix proxy) ─────────
# See header note: the true per-stage MACS peak_called_in field didn't
# survive into this RDS, so presence is read off the peak x cell counts
# matrix instead -- a peak counts toward a stage if accessible (count > 0)
# in at least MIN_CELLS_DETECTED cells of that stage.
MIN_CELLS_DETECTED <- 3
counts_mat <- GetAssayData(uro, assay = "ATAC", layer = "counts")

stage_summary <- lapply(STAGE_ORDER, function(st) {
  cells_st  <- colnames(uro)[uro$Age == st]
  detected  <- Matrix::rowSums(counts_mat[, cells_st, drop = FALSE] > 0) >= MIN_CELLS_DETECTED
  data.frame(Age = st, Category = peak_category[detected])
}) %>%
  bind_rows() %>%
  count(Age, Category, name = "n_peaks") %>%
  mutate(Age      = factor(Age, levels = STAGE_ORDER),
         Category = factor(Category, levels = ANNO_LEVELS)) %>%
  arrange(Age, Category)

stage_totals <- stage_summary %>%
  group_by(Age) %>%
  summarise(n_total = sum(n_peaks), .groups = "drop")

message(sprintf("  Peak detection: count > 0 in >= %d cells of a stage", MIN_CELLS_DETECTED))
message("  Peak annotation summary by stage:")
print(as.data.frame(stage_summary))
write.csv(stage_summary, file.path(OUT_DIR, "UrotheliumDevelopment_ATAC_PeakAnnotation_Summary.csv"),
          row.names = FALSE)

figI <- ggplot(stage_summary, aes(x = Age, y = n_peaks, fill = Category)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.3) +
  geom_text(data = stage_totals, aes(x = Age, y = n_total, label = n_total),
            inherit.aes = FALSE, vjust = -0.4, size = 3.2) +
  scale_fill_manual(values = ANNO_COLORS, name = "Genomic\nannotation") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Developmental stage", y = "Number of peaks",
       title = "Figure I. Peak-set genomic annotation by stage",
       subtitle = sprintf("Peaks detected per stage (count > 0 in >= %d cells); ChIPseeker vs. EnsDb.Mmusculus.v79",
                           MIN_CELLS_DETECTED)) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave(file.path(OUT_DIR, "FigI_UrotheliumDevelopment_ATAC_PeakAnnotation_Barplot.pdf"),
       figI, width = 7.5, height = 4.5)
message("  Saved: FigI_UrotheliumDevelopment_ATAC_PeakAnnotation_Barplot.pdf")

################################################################################
# Figure K: snATAC-seq TSS enrichment score distribution by stage
################################################################################
message("\n==> Building Figure K (TSS enrichment score violin by stage) ...")

tss_df <- data.frame(Age = uro$Age, TSS.enrichment = uro$TSS.enrichment)

message("  TSS.enrichment summary by stage:")
print(tss_df %>% group_by(Age) %>%
        summarise(median = median(TSS.enrichment), n = n(), .groups = "drop"))

figK <- ggplot(tss_df, aes(x = Age, y = TSS.enrichment, fill = Age)) +
  geom_violin(trim = TRUE, linewidth = 0.3, color = "black") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, linewidth = 0.3) +
  scale_fill_manual(values = STAGE_COLORS, guide = "none") +
  labs(x = "Developmental stage", y = "TSS enrichment score",
       title = "Figure K. snATAC-seq TSS enrichment score across development",
       subtitle = "Per-cell TSS.enrichment (Signac, fast mode); box = median/IQR") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave(file.path(OUT_DIR, "FigK_UrotheliumDevelopment_ATAC_TSSEnrichmentScore_ViolinByStage.pdf"),
       figK, width = 7.5, height = 4.5)
message("  Saved: FigK_UrotheliumDevelopment_ATAC_TSSEnrichmentScore_ViolinByStage.pdf")

message("\n==> Done.")
