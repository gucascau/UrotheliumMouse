################################################################################
# 04_UrotheliumDevelopment_CytoTRACE.R
#
# CytoTRACE differentiation score for the Uro-only population, as an
# independent (non-trajectory-graph-based) check on the slingshot ordering
# from 03_UrotheliumDevelopment_Pseudotime.R. CytoTRACE predicts relative
# differentiation status per cell from transcriptional diversity (roughly:
# less-differentiated cells express more genes and have more diffuse,
# covariant expression) rather than from a fitted low-dimensional curve, so
# agreement between the two is a useful cross-check.
#
# Raw counts:
#   The Chen2025 RDS object only ships a log-normalized "data" layer (see
#   [[project_chen2025_devatlas_structure]] memory note) -- no counts layer.
#   However the *source h5ad*
#   (UsedSingleCellsRawResults/RenalUrothelium/...Chen2025NatGenet.h5ad)
#   still carries the original integer UMI counts in its `raw.X` slot
#   (verified: same shape/gene order as `X`, integer-valued, max ~860).
#   04a_extract_uro_rawcounts_h5ad.py pulls just the 715 Uro-only rows out of
#   raw.X (no need to load/re-preprocess the full 203k-cell object) and
#   writes them to UrotheliumOnly_rawcounts.mtx. This script uses those true
#   raw counts if present; if 04a hasn't been run, it falls back to an
#   expm1()-recovered pseudo-count matrix from the log-normalized data (the
#   exact inverse of Seurat's log1p normalization) -- an approximation that
#   preserves each cell's gene-expression proportions and its zero/nonzero
#   detection pattern, but not its true absolute sequencing depth.
#
# Figure 7. Urothelium differentiation score (CytoTRACE)
#   (A) UMAP colored by differentiation score (1 - CytoTRACE stemness score,
#       so color direction matches the other figures: dark = immature,
#       bright = mature), with the slingshot curve overlaid.
#   (B) Differentiation score distribution by developmental stage.
#   (C) Differentiation score vs. slingshot pseudotime -- how well the two
#       independent methods agree on ordering (they do not agree strongly
#       when this was first run on the expm1-recovered pseudo-counts; see
#       whether true raw counts change that).
#
# Input:  UrotheliumDevelopmentScripts/output/UrotheliumOnly_pseudotime.rds
#         (written by 03_UrotheliumDevelopment_Pseudotime.R)
#         UrotheliumOnly_rawcounts.mtx (+ _genes.txt/_barcodes.txt), optional
#         (written by 04a_extract_uro_rawcounts_h5ad.py -- run it first for
#         true raw counts instead of the expm1 fallback)
# Output: Fig7_UrotheliumDevelopment_CytoTRACE.pdf
#         UrotheliumOnly_pseudotime_cytotrace.rds
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(CytoTRACE)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(viridisLite)
})

SCRIPT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/UrotheliumDevelopmentScripts"
OUT_DIR    <- file.path(SCRIPT_DIR, "output")

STAGE_ORDER  <- c("E16.5", "P0", "W3", "W12", "W52", "W92")
STAGE_COLORS <- setNames(viridisLite::viridis(length(STAGE_ORDER)), STAGE_ORDER)

# ── Load the pseudotime object ──────────────────────────────────────────────
message("==> Loading Uro-only object with pseudotime ...")
pt_rds <- file.path(OUT_DIR, "UrotheliumOnly_pseudotime.rds")
if (!file.exists(pt_rds)) {
  stop("Missing ", pt_rds, " -- run 03_UrotheliumDevelopment_Pseudotime.R first.")
}
uro <- readRDS(pt_rds)
uro$Age <- factor(uro$Age, levels = STAGE_ORDER)
umap_emb <- Embeddings(uro, "umap")

################################################################################
# CytoTRACE
################################################################################
raw_mtx <- file.path(OUT_DIR, "UrotheliumOnly_rawcounts.mtx")
raw_genes_file <- file.path(OUT_DIR, "UrotheliumOnly_rawcounts_genes.txt")
raw_barcodes_file <- file.path(OUT_DIR, "UrotheliumOnly_rawcounts_barcodes.txt")

if (all(file.exists(raw_mtx, raw_genes_file, raw_barcodes_file))) {
  message("\n==> Loading TRUE raw counts extracted from the source h5ad (raw.X) ...")
  pseudo_counts <- as.matrix(Matrix::readMM(raw_mtx))
  rownames(pseudo_counts) <- readLines(raw_genes_file)
  colnames(pseudo_counts) <- readLines(raw_barcodes_file)
  # Match to uro's cell order (04a writes them in the same order it read
  # them from UrotheliumOnly_barcodes.txt, which was uro's colnames -- this
  # match() is just a safety check, not expected to reorder anything).
  pseudo_counts <- pseudo_counts[, colnames(uro), drop = FALSE]
  message(sprintf("  %d genes x %d cells, value range [%.0f, %.0f], all-integer = %s",
                  nrow(pseudo_counts), ncol(pseudo_counts),
                  min(pseudo_counts), max(pseudo_counts),
                  isTRUE(all.equal(pseudo_counts, round(pseudo_counts)))))
} else {
  message("\n==> UrotheliumOnly_rawcounts.mtx not found -- run ",
          "04a_extract_uro_rawcounts_h5ad.py first for true raw counts. ",
          "Falling back to an expm1()-recovered pseudo-count matrix ",
          "(exact inverse of the log1p-normalized data layer; preserves ",
          "expression proportions and zero pattern, not true sequencing depth).")
  pseudo_counts <- as.matrix(expm1(GetAssayData(uro, assay = "originalexp", layer = "data")))
  message(sprintf("  %d genes x %d cells, value range [%.1f, %.1f]",
                  nrow(pseudo_counts), ncol(pseudo_counts),
                  min(pseudo_counts), max(pseudo_counts)))
}

message("\n==> Running CytoTRACE ...")
set.seed(1)
cyto_res <- CytoTRACE(mat = pseudo_counts, enableFast = FALSE)

if (length(cyto_res$filteredCells) > 0) {
  message(sprintf("  %d cell(s) filtered by CytoTRACE for low quality: %s",
                  length(cyto_res$filteredCells),
                  paste(head(cyto_res$filteredCells, 5), collapse = ", ")))
}

# Raw CytoTRACE score: 1 = least differentiated (stem-like), 0 = most
# differentiated. We also store 1 - score as "DifferentiationScore" so its
# color direction (dark = immature, bright = mature) matches every other
# figure in this series (Age, Pseudotime).
uro$CytoTRACE <- NA_real_
uro$CytoTRACE[names(cyto_res$CytoTRACE)] <- cyto_res$CytoTRACE
uro$DifferentiationScore <- 1 - uro$CytoTRACE

message("\n  Mean differentiation score by stage (sanity check -- should increase E16.5 -> W92):")
print(uro@meta.data %>% group_by(Age) %>%
        summarise(mean_score = mean(DifferentiationScore, na.rm = TRUE), .groups = "drop"))

rho_stage <- cor(uro$DifferentiationScore, as.integer(uro$Age),
                  method = "spearman", use = "complete.obs")
rho_pt    <- cor(uro$DifferentiationScore, uro$Pseudotime,
                  method = "spearman", use = "complete.obs")
message(sprintf("  Spearman correlation(differentiation score, chronological stage) = %.3f", rho_stage))
message(sprintf("  Spearman correlation(differentiation score, slingshot pseudotime) = %.3f", rho_pt))

# Top genes driving the CytoTRACE axis, mapped back to symbols for readability.
fm <- uro[["originalexp"]]@meta.features
top_cyto_genes <- head(cyto_res$cytoGenes, 15)
top_symbols <- fm$gene_symbols[match(names(top_cyto_genes), rownames(fm))]
message("\n  Top 15 genes positively correlated with CytoTRACE stemness score:")
print(data.frame(gene = top_symbols, correlation = round(top_cyto_genes, 3), row.names = NULL))

saveRDS(uro, file.path(OUT_DIR, "UrotheliumOnly_pseudotime_cytotrace.rds"))
message("\n  Saved: UrotheliumOnly_pseudotime_cytotrace.rds")

################################################################################
# Figure 7: differentiation score
################################################################################
message("\n==> Building Figure 7 (CytoTRACE differentiation score) ...")

# Re-derive the slingshot curve in UMAP space the same way Figure 5 did, so
# it can be overlaid here without re-fitting slingshot.
suppressPackageStartupMessages(library(slingshot))
pca_emb <- Embeddings(uro, "pca")
sds      <- slingshot(pca_emb, clusterLabels = uro$seurat_clusters, start.clus = "4")
sds_umap <- embedCurves(sds, umap_emb)
curve_df <- slingCurves(sds_umap, as.df = TRUE) %>% arrange(Lineage, Order)
arrow_spec <- arrow(length = unit(0.12, "inches"), type = "closed")

umap_df <- data.frame(umap_emb, DifferentiationScore = uro$DifferentiationScore)

p7a <- ggplot(umap_df, aes(x = umap_1, y = umap_2)) +
  geom_point(aes(color = DifferentiationScore), size = 1.3) +
  geom_path(data = curve_df, aes(x = umap_1, y = umap_2, group = Lineage),
            color = "black", linewidth = 0.8, arrow = arrow_spec) +
  scale_color_viridis_c(option = "plasma") +
  labs(color = "Differentiation\nscore") +
  ggtitle("CytoTRACE differentiation score") +
  coord_fixed() +
  theme_classic() +
  theme(aspect.ratio = 1, plot.title = element_text(face = "bold"),
        legend.title = element_text(size = 9))

p7b <- ggplot(uro@meta.data, aes(x = Age, y = DifferentiationScore, fill = Age)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
  scale_fill_manual(values = STAGE_COLORS, guide = "none") +
  labs(x = "Developmental stage", y = "Differentiation score") +
  ggtitle("By stage") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

p7c <- ggplot(uro@meta.data, aes(x = Pseudotime, y = DifferentiationScore)) +
  geom_point(aes(color = Age), size = 0.9, alpha = 0.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 4),
              color = "black", linewidth = 0.8, se = TRUE) +
  scale_color_manual(values = STAGE_COLORS) +
  labs(x = "Slingshot pseudotime", y = "Differentiation score", color = "Stage",
       subtitle = sprintf("Spearman rho = %.2f", rho_pt)) +
  ggtitle("vs. slingshot pseudotime") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

fig7 <- (p7a | p7b) / p7c +
  plot_annotation(title = "Figure 7. Urothelium differentiation score (CytoTRACE)",
                   theme = theme(plot.title = element_text(face = "bold", size = 14)))

ggsave(file.path(OUT_DIR, "Fig7_UrotheliumDevelopment_CytoTRACE.pdf"),
       fig7, width = 10, height = 10)
message("  Saved: Fig7_UrotheliumDevelopment_CytoTRACE.pdf")

message("\n==> Done.")
