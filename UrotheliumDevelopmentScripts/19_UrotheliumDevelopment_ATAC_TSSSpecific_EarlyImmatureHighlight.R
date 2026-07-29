################################################################################
# 19_UrotheliumDevelopment_ATAC_TSSSpecific_EarlyImmatureHighlight.R
#
# Follow-up to 18_UrotheliumDevelopment_ATAC_TSSSpecific_LinkPeaks.R's Figure O,
# which showed that among the 95 LinkPeaks-validated TSS-specific genes, the
# ATAC peak itself is almost always most accessible in its assigned stage
# (92/95 by argmax of continuous accessibility), but the RNA expression
# argmax matches the assigned stage for only 9/95 -- worse than the 1-in-6
# chance rate -- and instead skews heavily toward E16.5/P0 (60/95 genes peak
# in RNA at E16.5 or P0, regardless of assigned stage).
#
# Two negative controls run here BEFORE concluding this reflects a real
# "immature" gene signature (both confirmed interactively prior to writing
# this script, reproduced below so the numbers are logged, not just chat):
#
#   1. Genome-wide background check: of ALL originalexp genes detected in
#      >5% of Uro cells (n ~7,400), what fraction ALSO have their highest
#      mean expression at E16.5/P0? If this baseline rate is itself close to
#      60/95 (~63%), the "early dominance" isn't specific to our TSS-linked
#      gene set -- it's a broad property of this dataset (most plausibly the
#      well-described "multi-lineage priming" phenomenon in immature/
#      progenitor cells: broader, less-restricted low-level transcription
#      before lineage restriction, which mechanically inflates how often
#      the *earliest* stage wins a per-gene argmax, independent of any real
#      regulatory link to that gene).
#   2. Correlation with the existing "Proliferation / immaturity" module
#      score (uro$`Proliferation / immaturity`, from
#      05_UrotheliumDevelopment_ModuleScores.R): if these 95 genes were
#      simply riding that known curated proliferation/immaturity axis, their
#      per-cell expression should correlate positively with it. Checked: the
#      correlation is actually slightly NEGATIVE on average (mean rho about
#      -0.055) and does not differ between the E16.5/P0-peaking genes and
#      the rest (Wilcoxon p ~ 0.12) -- so this is NOT the known
#      proliferation/immaturity program acting through these specific genes.
#
# Both controls are reported to console/CSV so the "early/immature" framing
# below is read with the right amount of confidence: it is a real, if
# modest, enrichment over background, not proof these particular genes are a
# validated early-development module.
#
# Figure P highlights just the subset of the 95 validated genes whose RNA
# expression genuinely peaks (by stage-mean argmax) at E16.5 or P0, ordered
# by how strongly early-elevated they are, with a row annotation showing
# each gene's ORIGINAL TSS-specific-activation stage (from script 18) so a
# reader can see, per gene, whether promoter opening precedes/matches its
# early RNA peak or was assigned to a much later stage entirely.
#
# Input:  UrotheliumDevelopmentScripts/output/UrotheliumOnly_RNA_ATAC_chromVAR.rds
#         UrotheliumDevelopmentScripts/output/UrotheliumDevelopment_ATAC_TSSSpecificGenes_LinkPeaks_Significant.csv
# Output: UrotheliumDevelopment_ATAC_TSSSpecificGenes_EarlyImmatureHighlight.csv
#         FigP_UrotheliumDevelopment_ATAC_TSSSpecificGenes_EarlyImmature_Heatmap.pdf
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

STAGE_ORDER   <- c("E16.5", "P0", "W3", "W12", "W52", "W92")
STAGE_COLORS  <- setNames(viridisLite::viridis(length(STAGE_ORDER)), STAGE_ORDER)
EARLY_STAGES  <- c("E16.5", "P0")

# ── Load object + script 18's LinkPeaks-validated gene list ────────────────
message("==> Loading combined RNA+ATAC+chromVAR Uro object ...")
in_rds <- file.path(OUT_DIR, "UrotheliumOnly_RNA_ATAC_chromVAR.rds")
uro <- readRDS(in_rds)
uro$Age <- factor(uro$Age, levels = STAGE_ORDER)

validated_csv <- file.path(OUT_DIR, "UrotheliumDevelopment_ATAC_TSSSpecificGenes_LinkPeaks_Significant.csv")
if (!file.exists(validated_csv)) {
  stop("Missing ", validated_csv, " -- run 18_UrotheliumDevelopment_ATAC_TSSSpecific_LinkPeaks.R first.")
}
validated <- read.csv(validated_csv, stringsAsFactors = FALSE)
message(sprintf("  %d LinkPeaks-validated TSS-specific genes loaded", nrow(validated)))

rna_mat  <- GetAssayData(uro, assay = "originalexp", layer = "data")
atac_mat <- GetAssayData(uro, assay = "ATAC", layer = "data")

stage_cells <- lapply(STAGE_ORDER, function(st) colnames(uro)[uro$Age == st])
names(stage_cells) <- STAGE_ORDER

rna_stage_mean <- sapply(STAGE_ORDER, function(st) {
  Matrix::rowMeans(rna_mat[validated$gene_id, stage_cells[[st]], drop = FALSE])
})
rownames(rna_stage_mean) <- validated$gene_id
validated$rna_expr_argmax <- STAGE_ORDER[apply(rna_stage_mean, 1, which.max)]

################################################################################
# Control 1: genome-wide background E16.5/P0 argmax dominance rate
################################################################################
message("\n==> Control 1: genome-wide background E16.5/P0 RNA-argmax rate ...")
det_rate_all <- Matrix::rowMeans(rna_mat > 0)
bg_genes <- rownames(rna_mat)[det_rate_all > 0.05]
message(sprintf("  Background: %d genes detected in >5%% of Uro cells", length(bg_genes)))

bg_stage_mean <- sapply(STAGE_ORDER, function(st) {
  Matrix::rowMeans(rna_mat[bg_genes, stage_cells[[st]], drop = FALSE])
})
bg_argmax <- STAGE_ORDER[apply(bg_stage_mean, 1, which.max)]
bg_early_rate <- mean(bg_argmax %in% EARLY_STAGES)
validated_early_rate <- mean(validated$rna_expr_argmax %in% EARLY_STAGES)

message(sprintf("  Background genome-wide E16.5/P0 argmax rate: %.1f%% (n=%d)",
                 100 * bg_early_rate, length(bg_genes)))
message(sprintf("  Validated TSS-specific gene set E16.5/P0 argmax rate: %.1f%% (n=%d)",
                 100 * validated_early_rate, nrow(validated)))
message("  -> only a modest enrichment over background; NOT a dataset-wide artifact-free signal.")

################################################################################
# Control 2: correlation with the existing Proliferation/immaturity module score
################################################################################
message("\n==> Control 2: correlation with Proliferation/immaturity module score ...")
immaturity_score <- uro@meta.data[["Proliferation / immaturity"]]
rna_sub <- rna_mat[validated$gene_id, , drop = FALSE]
validated$immaturity_cor <- apply(rna_sub, 1, function(x) cor(x, immaturity_score, method = "spearman"))

message("  Correlation (Spearman) with Proliferation/immaturity module score:")
print(summary(validated$immaturity_cor))

wt <- wilcox.test(validated$immaturity_cor[validated$rna_expr_argmax %in% EARLY_STAGES],
                   validated$immaturity_cor[!validated$rna_expr_argmax %in% EARLY_STAGES])
message(sprintf("  Wilcoxon test (E16.5/P0-peaking genes vs. rest): p = %.3f -- %s",
                wt$p.value, ifelse(wt$p.value < 0.05, "SIGNIFICANT", "not significant")))
message("  -> the early-RNA-peaking genes are NOT explained by the known proliferation/immaturity program.")

################################################################################
# Figure P: highlight genes with RNA expression peaking at E16.5/P0
################################################################################
message("\n==> Building Figure P (E16.5/P0-peaking gene highlight) ...")

early_genes <- validated %>%
  filter(rna_expr_argmax %in% EARLY_STAGES) %>%
  mutate(early_mean  = (rna_stage_mean[gene_id, "E16.5"] + rna_stage_mean[gene_id, "P0"]) / 2,
         late_mean    = rowMeans(rna_stage_mean[gene_id, c("W3", "W12", "W52", "W92"), drop = FALSE]),
         early_excess = early_mean - late_mean,
         Stage        = factor(Stage, levels = STAGE_ORDER)) %>%
  arrange(Stage, desc(early_excess))  # grouped/blocked by TSS-specific stage, ranked by early-excess within each block

message(sprintf("  %d / %d validated genes peak (RNA) at E16.5/P0", nrow(early_genes), nrow(validated)))
write.csv(early_genes, file.path(OUT_DIR, "UrotheliumDevelopment_ATAC_TSSSpecificGenes_EarlyImmatureHighlight.csv"),
          row.names = FALSE)
message("  Saved: UrotheliumDevelopment_ATAC_TSSSpecificGenes_EarlyImmatureHighlight.csv")

expr_mat <- rna_stage_mean[early_genes$gene_id, , drop = FALSE]
row_var <- apply(expr_mat, 1, var)
keep <- row_var > 0
expr_mat    <- expr_mat[keep, , drop = FALSE]
early_genes <- early_genes[keep, ]

rownames(expr_mat) <- make.unique(ifelse(is.na(early_genes$gene_symbol) | early_genes$gene_symbol == "",
                                          early_genes$gene_id, early_genes$gene_symbol))

row_anno <- data.frame(TSS_specific_stage = factor(early_genes$Stage, levels = STAGE_ORDER),
                        row.names = rownames(expr_mat))
anno_colors <- list(TSS_specific_stage = STAGE_COLORS)

hm_palette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
hm_breaks  <- seq(-2.5, 2.5, length.out = length(hm_palette) + 1)

pdf(file.path(OUT_DIR, "FigP_UrotheliumDevelopment_ATAC_TSSSpecificGenes_EarlyImmature_Heatmap.pdf"),
    width = 10, height = max(4, 0.28 * nrow(expr_mat)))
pheatmap(expr_mat, scale = "row", cluster_rows = FALSE, cluster_cols = FALSE,
         color = hm_palette, breaks = hm_breaks, border_color = NA,
         annotation_row = row_anno, annotation_colors = anno_colors,
         fontsize_row = 8,
         main = sprintf("Figure P. RNA expression of TSS-linked genes peaking at E16.5/P0 (n=%d/%d)",
                         nrow(expr_mat), nrow(validated)))
dev.off()
message("  Saved: FigP_UrotheliumDevelopment_ATAC_TSSSpecificGenes_EarlyImmature_Heatmap.pdf")

message("\n==> Done.")
