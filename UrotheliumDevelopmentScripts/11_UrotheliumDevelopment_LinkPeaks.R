################################################################################
# 11_UrotheliumDevelopment_LinkPeaks.R
#
# Peak-gene cis-regulatory links (Signac::LinkPeaks) for Spp1 and a panel of
# barrier/differentiation markers, three figures:
#   Panel D: genome-browser coverage/track plot at the Spp1 locus across the
#     6 stages, linked peak(s) annotated, cross-referenced against Runx1
#     motif presence (from 10_UrotheliumDevelopment_chromVAR.R's motif scan).
#   Panel E: LinkPeaks score/significance summary across an 11-gene panel
#     (Spp1, Upk1a/1b/3a, Krt20, plus "novel" markers S100a6, Sprr1a, Psca,
#     Gsdmc2, Fxyd3, Snx31) -- are these genes under real cis-regulatory
#     control, not just incidentally co-expressed.
#   Panel F: stage-resolved accessibility (E16.5->W92) at the specific
#     Spp1-linked peak(s) -- the direct chromatin-accessibility parallel to
#     Spp1's RNA-expression-by-stage story elsewhere in this pipeline
#     (05's module scores, CellCommunicationScripts/02's CellChat run).
#
# genes.use restricts LinkPeaks to just these 11 genes -- this only limits
# which genes get correlation-tested, not the background peak universe
# LinkPeaks samples from for its null model (confirmed against Signac's
# behavior while planning this pipeline), so it's a cheap, statistically
# clean way to control cost, unlike restricting the underlying peak set.
#
# RegionStats (GC content / width per peak) was already computed once in
# script 08 and is required by LinkPeaks' background model -- not repeated
# here.
#
# Input:  UrotheliumDevelopmentScripts/output/UrotheliumOnly_RNA_ATAC_chromVAR.rds
#         (written by 10_UrotheliumDevelopment_chromVAR.R)
# Output: UrotheliumDevelopment_LinkPeaks_all.csv
#         UrotheliumDevelopment_LinkPeaks_GenePanel.csv
#         FigD_UrotheliumDevelopment_Spp1_CoveragePlot.pdf
#         FigE_UrotheliumDevelopment_LinkPeaks_GenePanel_Summary.pdf
#         FigF_UrotheliumDevelopment_Spp1LinkedPeak_Accessibility_ByStage.pdf
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(ggplot2)
  library(patchwork)
  library(viridisLite)
  library(dplyr)
})

SCRIPT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/UrotheliumDevelopmentScripts"
OUT_DIR    <- file.path(SCRIPT_DIR, "output")

STAGE_ORDER  <- c("E16.5", "P0", "W3", "W12", "W52", "W92")
STAGE_COLORS <- setNames(viridisLite::viridis(length(STAGE_ORDER)), STAGE_ORDER)

GENE_PANEL <- c("Spp1", "Upk1a", "Upk1b", "Upk3a", "Krt20",
                 "S100a6", "Sprr1a", "Psca", "Gsdmc2", "Fxyd3", "Snx31")

# ── Load ─────────────────────────────────────────────────────────────────────
message("==> Loading combined RNA+ATAC+chromVAR Uro object ...")
in_rds <- file.path(OUT_DIR, "UrotheliumOnly_RNA_ATAC_chromVAR.rds")
if (!file.exists(in_rds)) {
  stop("Missing ", in_rds, " -- run 10_UrotheliumDevelopment_chromVAR.R first.")
}
uro <- readRDS(in_rds)
uro$Age <- factor(uro$Age, levels = STAGE_ORDER)
DefaultAssay(uro) <- "ATAC"

fm_rna <- uro[["originalexp"]]@meta.features
missing_genes <- setdiff(GENE_PANEL, fm_rna$gene_symbols)
if (length(missing_genes) > 0) {
  stop("Gene(s) not found in originalexp gene_symbols: ", paste(missing_genes, collapse = ", "))
}

################################################################################
# LinkPeaks
################################################################################
message("\n==> Running LinkPeaks (genes.use = ", length(GENE_PANEL), " genes) ...")
uro <- LinkPeaks(
  object            = uro,
  peak.assay        = "ATAC",
  expression.assay  = "originalexp",
  genes.use         = GENE_PANEL
)

links <- Links(uro[["ATAC"]])
message(sprintf("  %d total peak-gene links found across the %d-gene panel", length(links), length(GENE_PANEL)))

links_df <- as.data.frame(links)
write.csv(links_df, file.path(OUT_DIR, "UrotheliumDevelopment_LinkPeaks_all.csv"), row.names = FALSE)
message("  Saved: UrotheliumDevelopment_LinkPeaks_all.csv")

sig_links <- links_df %>% filter(pvalue < 0.05) %>% arrange(gene, pvalue)
write.csv(sig_links, file.path(OUT_DIR, "UrotheliumDevelopment_LinkPeaks_GenePanel.csv"), row.names = FALSE)
message(sprintf("  %d significant links (p < 0.05):", nrow(sig_links)))
print(table(sig_links$gene))

################################################################################
# Panel D: Spp1 CoveragePlot with links, Runx1-motif cross-reference
################################################################################
message("\n==> Building Panel D (Spp1 CoveragePlot) ...")

spp1_links <- links_df %>% filter(gene == "Spp1")
region_highlight <- NULL
if (nrow(spp1_links) > 0) {
  spp1_peak_gr <- GRanges(spp1_links$seqnames, IRanges(spp1_links$start, spp1_links$end))
  # Cross-reference against Runx1 motif matches from script 10's AddMotifs scan.
  motif_matrix <- Motifs(uro[["ATAC"]])@data  # peaks x motifs, sparse logical/count
  motif_names  <- uro[["ATAC"]]@motifs@motif.names
  runx1_motif_id <- names(motif_names)[sapply(motif_names, function(nm) grepl("(^|::)Runx1(::|$)", nm, ignore.case = TRUE))]
  if (length(runx1_motif_id) > 0) {
    peak_names <- GRangesToString(granges(uro[["ATAC"]]))
    spp1_peak_names <- GRangesToString(spp1_peak_gr)
    has_runx1 <- spp1_peak_names %in% rownames(motif_matrix)[motif_matrix[, runx1_motif_id[1]] > 0]
    message(sprintf("  Spp1-linked peak(s) with a Runx1 motif: %d / %d", sum(has_runx1), length(has_runx1)))
    if (any(has_runx1)) {
      region_highlight <- spp1_peak_gr[has_runx1]
      write.csv(data.frame(peak = spp1_peak_names, has_Runx1_motif = has_runx1),
                file.path(OUT_DIR, "UrotheliumDevelopment_Spp1LinkedPeaks_Runx1Motif.csv"), row.names = FALSE)
    }
  } else {
    message("  Runx1 motif not found in the scanned motif set (check 10_UrotheliumDevelopment_chromVAR.R's TF coverage table).")
  }
} else {
  message("  No significant Spp1 links found -- CoveragePlot will still show accessibility/links track, just empty of Spp1-specific links. Reporting honestly rather than forcing a result.")
}

cov_args <- list(object = uro, region = "Spp1", features = "Spp1", assay = "ATAC",
                  expression.assay = "originalexp", group.by = "Age", links = TRUE,
                  extend.upstream = 5000, extend.downstream = 5000)
if (!is.null(region_highlight)) cov_args$region.highlight <- region_highlight
figD <- do.call(CoveragePlot, cov_args) &
  scale_fill_manual(values = STAGE_COLORS) &
  patchwork::plot_annotation(title = "Panel D. Spp1 locus: accessibility, expression, and linked peaks by stage")

ggsave(file.path(OUT_DIR, "FigD_UrotheliumDevelopment_Spp1_CoveragePlot.pdf"),
       figD, width = 9, height = 8)
message("  Saved: FigD_UrotheliumDevelopment_Spp1_CoveragePlot.pdf")

################################################################################
# Panel E: LinkPeaks summary across the gene panel
################################################################################
message("\n==> Building Panel E (LinkPeaks gene panel summary) ...")

gene_summary <- links_df %>%
  group_by(gene) %>%
  slice_min(pvalue, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  right_join(data.frame(gene = GENE_PANEL), by = "gene") %>%
  mutate(gene = factor(gene, levels = rev(GENE_PANEL)),
         neg_log10_p = -log10(pmax(pvalue, 1e-300)),
         significant = !is.na(pvalue) & pvalue < 0.05)

figE <- ggplot(gene_summary, aes(x = gene, y = score, fill = significant)) +
  geom_col(na.rm = TRUE) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "grey70"),
                     na.value = "grey90", labels = c("TRUE" = "p < 0.05", "FALSE" = "n.s."),
                     name = NULL) +
  labs(x = NULL, y = "Best peak-gene link score (Pearson, Signac LinkPeaks)",
       title = "Panel E. Peak-gene link strength across the marker panel",
       subtitle = "Best (lowest-p) link per gene; grey bars/blanks = no significant link found") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13))

ggsave(file.path(OUT_DIR, "FigE_UrotheliumDevelopment_LinkPeaks_GenePanel_Summary.pdf"),
       figE, width = 7, height = 5)
message("  Saved: FigE_UrotheliumDevelopment_LinkPeaks_GenePanel_Summary.pdf")

################################################################################
# Panel F: Spp1-linked peak accessibility by stage
################################################################################
message("\n==> Building Panel F (Spp1-linked peak accessibility by stage) ...")

if (nrow(spp1_links) == 0) {
  message("  No significant Spp1 link -- skipping Panel F (nothing to plot honestly).")
} else {
  best_spp1_peak <- spp1_links %>% slice_min(pvalue, n = 1, with_ties = FALSE) %>% pull(peak)
  acc <- GetAssayData(uro, assay = "ATAC", layer = "data")[best_spp1_peak, ]
  acc_df <- data.frame(Age = uro$Age, Accessibility = acc)

  stage_rank <- as.integer(uro$Age)
  rho_acc <- cor(acc, stage_rank, method = "spearman")

  figF <- ggplot(acc_df, aes(x = Age, y = Accessibility, fill = Age)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.7) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4,
             label = sprintf("rho = %.2f", rho_acc), size = 4) +
    scale_fill_manual(values = STAGE_COLORS, guide = "none") +
    labs(x = "Developmental stage", y = "Normalized ATAC accessibility",
         title = "Panel F. Spp1-linked peak accessibility across development",
         subtitle = sprintf("Peak: %s", best_spp1_peak)) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13))

  ggsave(file.path(OUT_DIR, "FigF_UrotheliumDevelopment_Spp1LinkedPeak_Accessibility_ByStage.pdf"),
         figF, width = 7, height = 6)
  message("  Saved: FigF_UrotheliumDevelopment_Spp1LinkedPeak_Accessibility_ByStage.pdf")
}

message("\n==> Done.")
