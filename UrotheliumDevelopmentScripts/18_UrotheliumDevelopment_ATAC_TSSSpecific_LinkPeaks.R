################################################################################
# 18_UrotheliumDevelopment_ATAC_TSSSpecific_LinkPeaks.R
#
# Figure O. RNA expression of TSS-specifically-activated genes, restricted to
#   genes whose promoter peak also has a statistically significant
#   Signac::LinkPeaks peak-gene correlation -- i.e. real evidence the peak's
#   accessibility tracks that specific gene's own expression across cells,
#   not just "peak sits within 3kb of this gene's TSS" (script 17's
#   nearest-TSS criterion) and not just "expression happens to vary by
#   stage" (no relationship to the peak required).
#
# Script 11 already ran LinkPeaks, but ONLY for an 11-gene curated urothelial
# marker panel (Spp1/Upk1a/Upk1b/Upk3a/Krt20/S100a6/Sprr1a/Psca/Gsdmc2/
# Fxyd3/Snx31) -- genes.use there restricted testing to exactly those genes.
# Checked interactively: the 231 TSS-specific genes from script 17
# (UrotheliumDevelopment_ATAC_TSSSpecificGenes_AllByStage.csv) share ZERO
# overlap with that 11-gene panel, so "restrict to genes with a significant
# script-11 LinkPeaks link" as literally asked would leave nothing to plot.
# This script instead reruns LinkPeaks with genes.use = the 230 unique
# TSS-specific gene symbols themselves (same method/defaults as script 11,
# just pointed at a different gene set) so every candidate actually gets a
# real peak-gene correlation test.
#
# Same recurring gotcha as script 11: LinkPeaks matches genes.use against
# rownames(GetAssayData(object, assay = expression.assay)), and Annotation()
# gene.coords lookups use gene SYMBOLS, not the Ensembl IDs that are
# originalexp's native rownames -- a small symbol-keyed assay is built here
# exactly as script 11 did, restricted to the TSS-specific gene set instead
# of the marker panel.
#
# Input:  UrotheliumDevelopmentScripts/output/UrotheliumOnly_RNA_ATAC_chromVAR.rds
#           (same input script 11 used -- already carries RegionStats +
#           Annotation, both required by LinkPeaks' background model)
#         UrotheliumDevelopmentScripts/output/UrotheliumDevelopment_ATAC_TSSSpecificGenes_AllByStage.csv
#           (written by 17_UrotheliumDevelopment_ATAC_TSSSpecific_Expression.R)
# Output: UrotheliumDevelopment_ATAC_TSSSpecificGenes_LinkPeaks_all.csv
#         UrotheliumDevelopment_ATAC_TSSSpecificGenes_LinkPeaks_Significant.csv
#         FigO_UrotheliumDevelopment_ATAC_TSSSpecificGenes_LinkPeaks_Expression_Heatmap.pdf
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
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

# ── Load the combined RNA+ATAC+chromVAR Uro object (same as script 11) ─────
message("==> Loading combined RNA+ATAC+chromVAR Uro object ...")
in_rds <- file.path(OUT_DIR, "UrotheliumOnly_RNA_ATAC_chromVAR.rds")
if (!file.exists(in_rds)) {
  stop("Missing ", in_rds, " -- run 10_UrotheliumDevelopment_chromVAR.R first.")
}
uro <- readRDS(in_rds)
uro$Age <- factor(uro$Age, levels = STAGE_ORDER)
DefaultAssay(uro) <- "ATAC"

# ── Load script 17's TSS-specific gene list ─────────────────────────────────
tss_csv <- file.path(OUT_DIR, "UrotheliumDevelopment_ATAC_TSSSpecificGenes_AllByStage.csv")
if (!file.exists(tss_csv)) {
  stop("Missing ", tss_csv, " -- run 17_UrotheliumDevelopment_ATAC_TSSSpecific_Expression.R first.")
}
tss_specific <- read.csv(tss_csv, stringsAsFactors = FALSE) %>%
  filter(in_rna_assay, !is.na(gene_symbol), gene_symbol != "")
message(sprintf("  %d stage x gene rows, %d unique gene symbols to test",
                 nrow(tss_specific), length(unique(tss_specific$gene_symbol))))

fm_rna <- uro[["originalexp"]]@meta.features
gene_symbols_use <- unique(tss_specific$gene_symbol)
missing_genes <- setdiff(gene_symbols_use, fm_rna$gene_symbols)
if (length(missing_genes) > 0) {
  message(sprintf("  Dropping %d symbol(s) not resolvable in gene_symbols (dup/ambiguous): %s",
                   length(missing_genes), paste(missing_genes, collapse = ", ")))
  gene_symbols_use <- setdiff(gene_symbols_use, missing_genes)
}

################################################################################
# LinkPeaks, scoped to the TSS-specific candidate gene set
################################################################################
message(sprintf("\n==> Running LinkPeaks (genes.use = %d TSS-specific genes) ...", length(gene_symbols_use)))

gene_ensembl_ids <- setNames(rownames(fm_rna)[match(gene_symbols_use, fm_rna$gene_symbols)], gene_symbols_use)
rna_panel_data <- GetAssayData(uro, assay = "originalexp", layer = "data")[gene_ensembl_ids, ]
rownames(rna_panel_data) <- names(gene_ensembl_ids)
uro[["rna_tss_specific"]] <- CreateAssayObject(data = rna_panel_data)

uro <- LinkPeaks(
  object            = uro,
  peak.assay        = "ATAC",
  expression.assay  = "rna_tss_specific",
  genes.use         = gene_symbols_use
)

links <- Links(uro[["ATAC"]])
message(sprintf("  %d total peak-gene links found across the %d TSS-specific genes", length(links), length(gene_symbols_use)))

links_df <- as.data.frame(links)
write.csv(links_df, file.path(OUT_DIR, "UrotheliumDevelopment_ATAC_TSSSpecificGenes_LinkPeaks_all.csv"),
          row.names = FALSE)
message("  Saved: UrotheliumDevelopment_ATAC_TSSSpecificGenes_LinkPeaks_all.csv")

sig_links <- links_df %>% filter(pvalue < 0.05) %>% arrange(gene, pvalue)
sig_genes <- unique(sig_links$gene)
message(sprintf("  %d significant peak-gene links (p < 0.05), covering %d / %d genes",
                 nrow(sig_links), length(sig_genes), length(gene_symbols_use)))

# ── Merge back onto the TSS-specific stage assignment, keep only significant genes ─
validated <- tss_specific %>%
  filter(gene_symbol %in% sig_genes) %>%
  left_join(sig_links %>% group_by(gene) %>% slice_min(pvalue, n = 1, with_ties = FALSE) %>%
              ungroup() %>% select(gene, link_score = score, link_pvalue = pvalue),
            by = c("gene_symbol" = "gene")) %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%
  mutate(Stage = factor(Stage, levels = STAGE_ORDER)) %>%
  arrange(Stage, link_pvalue)

message(sprintf("\n  Stage-specific TSS-activation genes WITH a significant LinkPeaks correlation: %d", nrow(validated)))
print(validated %>% count(Stage))
write.csv(validated, file.path(OUT_DIR, "UrotheliumDevelopment_ATAC_TSSSpecificGenes_LinkPeaks_Significant.csv"),
          row.names = FALSE)
message("  Saved: UrotheliumDevelopment_ATAC_TSSSpecificGenes_LinkPeaks_Significant.csv")

################################################################################
# Figure O: RNA expression of LinkPeaks-validated TSS-specific genes
################################################################################
if (nrow(validated) == 0) {
  message("\n  No genes survived the LinkPeaks significance filter -- skipping Figure O (nothing to plot honestly).")
} else {
  message("\n==> Building Figure O (expression heatmap, LinkPeaks-validated genes) ...")

  data_mat <- GetAssayData(uro, assay = "originalexp", layer = "data")
  expr_mat <- sapply(STAGE_ORDER, function(st) {
    cells_st <- colnames(uro)[uro$Age == st]
    Matrix::rowMeans(data_mat[validated$gene_id, cells_st, drop = FALSE])
  })
  rownames(expr_mat) <- validated$gene_id
  colnames(expr_mat) <- STAGE_ORDER

  row_var <- apply(expr_mat, 1, var)
  keep <- row_var > 0
  message(sprintf("  Dropping %d gene(s) with zero expression variance across stages", sum(!keep)))
  expr_mat  <- expr_mat[keep, , drop = FALSE]
  validated <- validated[keep, ]

  rownames(expr_mat) <- make.unique(validated$gene_symbol)

  row_anno <- data.frame(TSS_specific_stage = validated$Stage, row.names = rownames(expr_mat))
  anno_colors <- list(TSS_specific_stage = STAGE_COLORS)

  hm_palette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
  hm_breaks  <- seq(-2.5, 2.5, length.out = length(hm_palette) + 1)

  pdf(file.path(OUT_DIR, "FigO_UrotheliumDevelopment_ATAC_TSSSpecificGenes_LinkPeaks_Expression_Heatmap.pdf"),
      width = 10, height = max(4, 0.28 * nrow(expr_mat)))
  pheatmap(expr_mat, scale = "row", cluster_rows = FALSE, cluster_cols = FALSE,
           color = hm_palette, breaks = hm_breaks, border_color = NA,
           annotation_row = row_anno, annotation_colors = anno_colors,
           fontsize_row = 8,
           main = sprintf("Figure O. RNA expression of TSS-specific + LinkPeaks-validated genes (n=%d)",
                           nrow(expr_mat)))
  dev.off()
  message("  Saved: FigO_UrotheliumDevelopment_ATAC_TSSSpecificGenes_LinkPeaks_Expression_Heatmap.pdf")
}

message("\n==> Done.")
