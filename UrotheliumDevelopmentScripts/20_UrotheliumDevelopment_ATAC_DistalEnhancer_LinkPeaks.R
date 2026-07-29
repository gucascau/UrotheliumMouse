################################################################################
# 20_UrotheliumDevelopment_ATAC_DistalEnhancer_LinkPeaks.R
#
# Distal-peak analog of 17/18's promoter analysis: identify stage-specific
# DISTAL/intergenic peaks (the largest annotation category, 37.8% of all
# peaks -- Fig I) and test whether any are under real cis-regulatory control
# of a nearby gene, using Signac::LinkPeaks rather than nearest-gene
# proximity to make the peak-to-gene assignment.
#
# Distal peaks cannot be assigned to "their" gene the way promoter peaks
# were (nearest TSS within 3kb) -- a distal element can regulate a gene tens
# to hundreds of kb away, skipping closer genes entirely. So instead of
# ChIPseeker's nearest-gene call, LinkPeaks itself (which already tests every
# peak within its search window against each candidate gene's expression) is
# used as the peak-to-gene assignment mechanism: run LinkPeaks broadly, then
# keep only the resulting significant links whose PEAK happens to be one of
# our stage-specific distal peaks.
#
# Candidate gene pool: genes within LinkPeaks' own default 500kb search
# window of any stage-specific distal peak (checked interactively: this is
# 2,299 well-detected genes, vs. 7,381 genome-wide -- a real but modest
# reduction, since a 500kb window already covers roughly a third of the
# well-detected genome). Restricting to this window is not a shortcut that
# loses anything LinkPeaks could have found anyway -- it's exactly what
# LinkPeaks would search, just computed once up front to size genes.use.
#
# Same peak-presence proxy as scripts 14/15/17 (count > 0 in >= 3 cells of a
# stage, degree == 1 = stage-specific), same ChIPseeker Distal/Downstream
# collapse used in Figure I, same symbol-keyed-assay pattern for LinkPeaks
# as scripts 11/18 (genes.use is matched against expression.assay rownames,
# which must be gene symbols to match Annotation()'s gene_name field).
#
# Input:  UrotheliumDevelopmentScripts/output/UrotheliumOnly_RNA_ATAC_chromVAR.rds
# Output: UrotheliumDevelopment_ATAC_DistalPeaks_LinkPeaks_all.csv
#         UrotheliumDevelopment_ATAC_DistalPeaks_LinkPeaks_Significant.csv
#         FigQ_UrotheliumDevelopment_ATAC_DistalEnhancerGenes_Expression_Heatmap.pdf
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(EnsDb.Mmusculus.v79)
  library(ChIPseeker)
  library(Matrix)
  library(dplyr)
  library(pheatmap)
  library(RColorBrewer)
  library(viridisLite)
})

SCRIPT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/UrotheliumDevelopmentScripts"
OUT_DIR    <- file.path(SCRIPT_DIR, "output")

STAGE_ORDER  <- c("E16.5", "P0", "W3", "W12", "W52", "W92")
STAGE_COLORS <- setNames(viridisLite::viridis(length(STAGE_ORDER)), STAGE_ORDER)
MIN_CELLS_DETECTED <- 3
LINKPEAKS_WINDOW    <- 500000  # Signac::LinkPeaks default `distance`

# ── Load combined RNA+ATAC+chromVAR Uro object (same as scripts 11/18) ─────
message("==> Loading combined RNA+ATAC+chromVAR Uro object ...")
uro <- readRDS(file.path(OUT_DIR, "UrotheliumOnly_RNA_ATAC_chromVAR.rds"))
uro$Age <- factor(uro$Age, levels = STAGE_ORDER)
DefaultAssay(uro) <- "ATAC"

# ── Peak annotation + presence matrix (same as scripts 14/17) ──────────────
message("\n==> Annotating peaks and identifying stage-specific distal peaks ...")
peaks <- granges(uro[["ATAC"]])
peaks_ens <- peaks
seqlevelsStyle(peaks_ens) <- "Ensembl"
anno <- annotatePeak(peaks_ens, TxDb = EnsDb.Mmusculus.v79,
                      tssRegion = c(-3000, 3000), verbose = FALSE)
anno_df <- as.data.frame(anno)
is_distal <- grepl("Distal|Downstream", anno_df$annotation)

counts_mat <- GetAssayData(uro, assay = "ATAC", layer = "counts")
presence_mat <- sapply(STAGE_ORDER, function(st) {
  cells_st <- colnames(uro)[uro$Age == st]
  Matrix::rowSums(counts_mat[, cells_st, drop = FALSE] > 0) >= MIN_CELLS_DETECTED
})
degree <- rowSums(presence_mat)

distal_specific <- is_distal & degree == 1
stage_of_peak <- STAGE_ORDER[apply(presence_mat, 1, function(x) if (any(x)) which(x)[1] else NA)]
message(sprintf("  %d stage-specific distal peaks (of %d distal, %d total)",
                 sum(distal_specific), sum(is_distal), length(peaks)))
print(table(stage_of_peak[distal_specific]))

distal_peak_names  <- rownames(anno_df)[distal_specific]
distal_peaks_ens    <- peaks_ens[distal_specific]

# ── Candidate genes.use: genes within LinkPeaks' own 500kb window ──────────
message(sprintf("\n==> Building candidate gene pool (genes within %dkb of a stage-specific distal peak) ...",
                 LINKPEAKS_WINDOW / 1000))
genes_gr <- genes(EnsDb.Mmusculus.v79)
hits <- findOverlaps(distal_peaks_ens + LINKPEAKS_WINDOW, genes_gr)
candidate_gene_ids <- unique(genes_gr$gene_id[subjectHits(hits)])

fm_rna <- uro[["originalexp"]]@meta.features
det_rate <- Matrix::rowMeans(GetAssayData(uro, assay = "originalexp", layer = "data") > 0)
well_detected_ids <- rownames(fm_rna)[det_rate > 0.05]
candidate_gene_ids <- intersect(candidate_gene_ids, well_detected_ids)
gene_symbols_use <- unique(na.omit(fm_rna$gene_symbols[match(candidate_gene_ids, rownames(fm_rna))]))
message(sprintf("  %d candidate genes (well-detected, resolvable symbol)", length(gene_symbols_use)))

################################################################################
# LinkPeaks, scoped to the distal-peak neighborhood gene pool
################################################################################
message(sprintf("\n==> Running LinkPeaks (genes.use = %d genes near stage-specific distal peaks) ...",
                 length(gene_symbols_use)))
gene_ensembl_ids <- setNames(rownames(fm_rna)[match(gene_symbols_use, fm_rna$gene_symbols)], gene_symbols_use)
rna_panel_data <- GetAssayData(uro, assay = "originalexp", layer = "data")[gene_ensembl_ids, ]
rownames(rna_panel_data) <- names(gene_ensembl_ids)
uro[["rna_distal_candidates"]] <- CreateAssayObject(data = rna_panel_data)

uro <- LinkPeaks(
  object            = uro,
  peak.assay        = "ATAC",
  expression.assay  = "rna_distal_candidates",
  genes.use         = gene_symbols_use,
  distance          = LINKPEAKS_WINDOW
)

links_df <- as.data.frame(Links(uro[["ATAC"]]))
message(sprintf("  %d total peak-gene links found across the %d candidate genes", nrow(links_df), length(gene_symbols_use)))
write.csv(links_df, file.path(OUT_DIR, "UrotheliumDevelopment_ATAC_DistalPeaks_LinkPeaks_all.csv"), row.names = FALSE)
message("  Saved: UrotheliumDevelopment_ATAC_DistalPeaks_LinkPeaks_all.csv")

# ── Keep only links whose PEAK is one of our stage-specific distal peaks ───
validated <- links_df %>%
  filter(pvalue < 0.05, peak %in% distal_peak_names) %>%
  mutate(Stage = stage_of_peak[match(peak, rownames(anno_df))]) %>%
  group_by(gene) %>%
  slice_min(pvalue, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  left_join(fm_rna %>% mutate(gene_id = rownames(fm_rna)) %>% select(gene_id, gene_symbols),
            by = c("gene" = "gene_symbols")) %>%
  rename(gene_symbol = gene, link_score = score, link_pvalue = pvalue) %>%
  mutate(Stage = factor(Stage, levels = STAGE_ORDER)) %>%
  arrange(Stage, link_pvalue)

message(sprintf("\n  Stage-specific DISTAL peaks with a significant LinkPeaks-assigned gene: %d", nrow(validated)))
print(validated %>% count(Stage))
write.csv(validated, file.path(OUT_DIR, "UrotheliumDevelopment_ATAC_DistalPeaks_LinkPeaks_Significant.csv"),
          row.names = FALSE)
message("  Saved: UrotheliumDevelopment_ATAC_DistalPeaks_LinkPeaks_Significant.csv")

################################################################################
# Figure Q: RNA expression of LinkPeaks-validated distal-enhancer target genes
################################################################################
if (nrow(validated) == 0) {
  message("\n  No genes survived -- skipping Figure Q (nothing to plot honestly).")
} else {
  message("\n==> Building Figure Q (expression heatmap, distal-enhancer-validated genes) ...")

  data_mat <- GetAssayData(uro, assay = "originalexp", layer = "data")
  expr_mat <- sapply(STAGE_ORDER, function(st) {
    cells_st <- colnames(uro)[uro$Age == st]
    Matrix::rowMeans(data_mat[validated$gene_id, cells_st, drop = FALSE])
  })
  rownames(expr_mat) <- validated$gene_id
  colnames(expr_mat) <- STAGE_ORDER

  rna_argmax <- STAGE_ORDER[apply(expr_mat, 1, which.max)]
  validated$rna_expr_argmax <- rna_argmax
  message(sprintf("  RNA expr argmax matches assigned (distal peak) stage: %d / %d",
                   sum(as.character(validated$Stage) == validated$rna_expr_argmax), nrow(validated)))
  write.csv(validated, file.path(OUT_DIR, "UrotheliumDevelopment_ATAC_DistalPeaks_LinkPeaks_Significant.csv"),
            row.names = FALSE)

  row_var <- apply(expr_mat, 1, var)
  keep <- row_var > 0
  message(sprintf("  Dropping %d gene(s) with zero expression variance across stages", sum(!keep)))
  expr_mat  <- expr_mat[keep, , drop = FALSE]
  validated <- validated[keep, ]

  rownames(expr_mat) <- make.unique(ifelse(is.na(validated$gene_symbol) | validated$gene_symbol == "",
                                            validated$gene_id, validated$gene_symbol))

  row_anno <- data.frame(distal_peak_stage = validated$Stage, row.names = rownames(expr_mat))
  anno_colors <- list(distal_peak_stage = STAGE_COLORS)

  hm_palette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
  hm_breaks  <- seq(-2.5, 2.5, length.out = length(hm_palette) + 1)

  pdf(file.path(OUT_DIR, "FigQ_UrotheliumDevelopment_ATAC_DistalEnhancerGenes_Expression_Heatmap.pdf"),
      width = 10, height = max(4, 0.28 * nrow(expr_mat)))
  pheatmap(expr_mat, scale = "row", cluster_rows = FALSE, cluster_cols = FALSE,
           color = hm_palette, breaks = hm_breaks, border_color = NA,
           annotation_row = row_anno, annotation_colors = anno_colors,
           fontsize_row = 8,
           main = sprintf("Figure Q. RNA expression of distal-enhancer-linked genes (n=%d)", nrow(expr_mat)))
  dev.off()
  message("  Saved: FigQ_UrotheliumDevelopment_ATAC_DistalEnhancerGenes_Expression_Heatmap.pdf")
}

message("\n==> Done.")
