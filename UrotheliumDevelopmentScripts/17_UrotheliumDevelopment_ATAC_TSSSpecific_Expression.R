################################################################################
# 17_UrotheliumDevelopment_ATAC_TSSSpecific_Expression.R
#
# Figure N. Stage-specific TSS/promoter activation -> RNA expression check.
#
# "TSS-specific activation" is operationalized here as: a ChIPseeker
# Promoter-annotated peak (same tssRegion = +/-3kb call as Figure I in
# 14_UrotheliumDevelopment_ATAC_QC_Figures.R) that is present in EXACTLY ONE
# stage (degree == 1 in the same presence/detection framework used for
# Figures I/J -- count > 0 in >= MIN_CELLS_DETECTED cells of that stage, and
# absent from every other stage by that same criterion). This reuses the
# exact peak-presence matrix logic from scripts 14/15 rather than
# introducing a new threshold.
#
# The Promoter peak's ChIPseeker geneId (Ensembl mouse ID, since
# EnsDb.Mmusculus.v79 was used as the TxDb) is taken as the "TSS-specifically
# activated gene" for that stage. geneId happens to match rownames(uro[["
# originalexp"]]) directly (both Ensembl IDs, confirmed interactively) --
# no symbol-based join needed to pull RNA expression for these genes.
#
# For the expression check, per stage the up-to-20 genes with the smallest
# |distanceToTSS| are kept for the heatmap (this pipeline's own convention
# for "top N per group" marker tables -- see
# UrotheliumDevelopment_StageMarkers_top20.csv from script 07); the full,
# unfiltered gene list (every stage-specific Promoter peak's gene) is saved
# separately with no cap.
#
# Expression is mean log-normalized "data" (originalexp assay, Uro cells
# only) per stage, z-scored per gene across the 6 stages -- NOT
# Seurat::AverageExpression(), which would exponentiate/back-transform the
# already-log data by default and misrepresent this as a linear-scale mean.
#
# Input:  UrotheliumDevelopmentScripts/output/UrotheliumOnly_RNA_ATAC.rds
# Output: UrotheliumDevelopment_ATAC_TSSSpecificGenes_AllByStage.csv (every
#           stage-specific promoter-peak gene, no cap)
#         UrotheliumDevelopment_ATAC_TSSSpecificGenes_Top20PerStage.csv
#           (subset plotted in Figure N)
#         FigN_UrotheliumDevelopment_ATAC_TSSSpecificGenes_Expression_Heatmap.pdf
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
MIN_CELLS_DETECTED <- 3  # same threshold as scripts 14/15
TOP_N_PER_STAGE    <- 20 # same "top N marker" convention as script 07

# ── Load the combined RNA+ATAC Uro object ───────────────────────────────────
message("==> Loading Uro RNA+ATAC object ...")
in_rds <- file.path(OUT_DIR, "UrotheliumOnly_RNA_ATAC.rds")
if (!file.exists(in_rds)) {
  stop("Missing ", in_rds, " -- run 08_UrotheliumDevelopment_ATAC_Build.R first.")
}
uro <- readRDS(in_rds)
uro$Age <- factor(uro$Age, levels = STAGE_ORDER)
message(sprintf("  %d cells, %d ATAC peaks", ncol(uro), nrow(uro[["ATAC"]])))

# ── Peak presence per stage (same proxy as Figures I/J) ────────────────────
message(sprintf("\n==> Peak presence matrix (count > 0 in >= %d cells/stage) ...", MIN_CELLS_DETECTED))
counts_mat <- GetAssayData(uro, assay = "ATAC", layer = "counts")
presence_mat <- sapply(STAGE_ORDER, function(st) {
  cells_st <- colnames(uro)[uro$Age == st]
  Matrix::rowSums(counts_mat[, cells_st, drop = FALSE] > 0) >= MIN_CELLS_DETECTED
})
rownames(presence_mat) <- rownames(counts_mat)
degree <- rowSums(presence_mat)

# ── Peak -> gene / genomic annotation (ChIPseeker vs. EnsDb.Mmusculus.v79) ──
message("\n==> Annotating peaks (ChIPseeker vs. EnsDb.Mmusculus.v79) ...")
peaks <- granges(uro[["ATAC"]])
peaks_ens <- peaks
seqlevelsStyle(peaks_ens) <- "Ensembl"
anno <- annotatePeak(peaks_ens, TxDb = EnsDb.Mmusculus.v79,
                      tssRegion = c(-3000, 3000), verbose = FALSE)
anno_df <- as.data.frame(anno)
is_promoter <- grepl("^Promoter", anno_df$annotation)
message(sprintf("  %d / %d peaks are Promoter-annotated", sum(is_promoter), length(peaks)))

# ── Stage-specific Promoter peaks -> genes ─────────────────────────────────
message("\n==> Identifying stage-specific TSS/promoter activation ...")
fm <- uro[["originalexp"]]@meta.features

tss_specific <- lapply(STAGE_ORDER, function(st) {
  idx <- which(is_promoter & degree == 1 & presence_mat[, st])
  if (length(idx) == 0) return(NULL)
  data.frame(
    Stage         = st,
    gene_id       = anno_df$geneId[idx],
    gene_symbol   = fm$gene_symbols[match(anno_df$geneId[idx], rownames(fm))],
    peak          = rownames(anno_df)[idx],
    distanceToTSS = anno_df$distanceToTSS[idx],
    stringsAsFactors = FALSE
  )
}) %>% bind_rows() %>%
  distinct(Stage, gene_id, .keep_all = TRUE) %>%   # a gene can have >1 stage-specific promoter peak
  mutate(in_rna_assay = gene_id %in% rownames(fm)) %>%
  arrange(factor(Stage, levels = STAGE_ORDER), abs(distanceToTSS))

message("  Stage-specific Promoter-peak genes (total vs. found in RNA assay):")
print(tss_specific %>% count(Stage, in_rna_assay) %>% tidyr::pivot_wider(names_from = in_rna_assay, values_from = n, values_fill = 0))

write.csv(tss_specific, file.path(OUT_DIR, "UrotheliumDevelopment_ATAC_TSSSpecificGenes_AllByStage.csv"),
          row.names = FALSE)
message("  Saved: UrotheliumDevelopment_ATAC_TSSSpecificGenes_AllByStage.csv")

# ── Top N per stage (closest to true TSS), restricted to genes in RNA assay ─
top_genes <- tss_specific %>%
  filter(in_rna_assay) %>%
  group_by(Stage) %>%
  slice_min(order_by = abs(distanceToTSS), n = TOP_N_PER_STAGE, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(Stage = factor(Stage, levels = STAGE_ORDER)) %>%
  arrange(Stage, abs(distanceToTSS))

message(sprintf("\n  Top-%d-per-stage genes selected for Figure N: %d genes total", TOP_N_PER_STAGE, nrow(top_genes)))
write.csv(top_genes, file.path(OUT_DIR, "UrotheliumDevelopment_ATAC_TSSSpecificGenes_Top20PerStage.csv"),
          row.names = FALSE)
message("  Saved: UrotheliumDevelopment_ATAC_TSSSpecificGenes_Top20PerStage.csv")

################################################################################
# Figure N: RNA expression of TSS-specifically-activated genes, by stage
################################################################################
message("\n==> Building Figure N (expression heatmap) ...")

data_mat <- GetAssayData(uro, assay = "originalexp", layer = "data")
expr_mat <- sapply(STAGE_ORDER, function(st) {
  cells_st <- colnames(uro)[uro$Age == st]
  Matrix::rowMeans(data_mat[top_genes$gene_id, cells_st, drop = FALSE])
})
rownames(expr_mat) <- top_genes$gene_id
colnames(expr_mat) <- STAGE_ORDER

# Drop genes with zero variance across stages (can't be z-scored, would be
# a flat/undefined row in the heatmap).
row_var <- apply(expr_mat, 1, var)
keep <- row_var > 0
message(sprintf("  Dropping %d gene(s) with zero expression variance across stages", sum(!keep)))
expr_mat  <- expr_mat[keep, , drop = FALSE]
top_genes <- top_genes[keep, ]

# Row labels: gene symbol, falling back to Ensembl ID if no symbol.
row_labels <- ifelse(is.na(top_genes$gene_symbol) | top_genes$gene_symbol == "",
                      top_genes$gene_id, top_genes$gene_symbol)
rownames(expr_mat) <- make.unique(row_labels)

row_anno <- data.frame(TSS_specific_stage = top_genes$Stage, row.names = rownames(expr_mat))
anno_colors <- list(TSS_specific_stage = STAGE_COLORS)

hm_palette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
hm_breaks  <- seq(-2.5, 2.5, length.out = length(hm_palette) + 1)

pdf(file.path(OUT_DIR, "FigN_UrotheliumDevelopment_ATAC_TSSSpecificGenes_Expression_Heatmap.pdf"),
    width = 10, height = max(4, 0.22 * nrow(expr_mat)))
pheatmap(expr_mat, scale = "row", cluster_rows = FALSE, cluster_cols = FALSE,
         color = hm_palette, breaks = hm_breaks, border_color = NA,
         annotation_row = row_anno, annotation_colors = anno_colors,
         fontsize_row = 6,
         main = sprintf("Figure N. RNA expression of TSS-specifically-activated genes (top %d/stage, n=%d)",
                         TOP_N_PER_STAGE, nrow(expr_mat)))
dev.off()
message("  Saved: FigN_UrotheliumDevelopment_ATAC_TSSSpecificGenes_Expression_Heatmap.pdf")

message("\n==> Done.")
