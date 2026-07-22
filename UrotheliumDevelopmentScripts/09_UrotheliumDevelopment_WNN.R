################################################################################
# 09_UrotheliumDevelopment_WNN.R
#
# Joint RNA+ATAC (WNN) analysis of the 715 Uro cells, two purposes:
#   Panel A: a WNN UMAP colored by stage -- sanity check that chromatin state
#     (or RNA+ATAC combined) separates cells by developmental stage similarly
#     to RNA alone (Fig1/Fig4/Fig5's stage-colored UMAPs).
#   Panel B: does an independent, ATAC-informed trajectory agree with the
#     existing RNA-only Slingshot trajectory (03_UrotheliumDevelopment_
#     Pseudotime.R)? This is the multiome answer to the CytoTRACE-vs-Slingshot
#     disagreement flagged in 04_UrotheliumDevelopment_CytoTRACE.R's header
#     (rho ~ 0.14 there) -- a THIRD, independent trajectory estimate that
#     either corroborates the RNA Slingshot direction or doesn't.
#
# WNN pseudotime design (mirrors 03's "fit on PCA, not 2D UMAP" rigor):
# FindMultiModalNeighbors doesn't produce a joint linear embedding (WNN is a
# weighted neighbor graph + per-cell modality weights, not a coordinate
# space) -- the only coordinate embedding it gives directly is a 2D UMAP.
# Fitting Slingshot straight on a 2D UMAP is a real step down in rigor (2D
# UMAP is known to introduce spurious loop/branch topology) versus this
# project's own RNA-only convention of fitting on ~15-dim PCA and only
# re-embedding the fitted curve into 2D for plotting. So here: a SEPARATE
# higher-dimensional WNN UMAP (10 components) is built purely for Slingshot
# fitting, distinct from the 2D wnn.umap used for Panel A / plotting -- same
# two-embeddings-for-two-purposes split 03 already uses (pca vs. umap).
#
# ATAC "lsi" reduction drops component 1 (RunTFIDF's first SVD component is
# near-universally a sequencing-depth/total-counts artifact, standard Signac
# practice, not a biological axis) -- dims.list starts LSI at 2.
# RNA "pca" reduction reused as-is from 02_UrotheliumDevelopment_Subcluster_
# UMAP.R (Uro-specific, 15 PCs) -- NOT the full-atlas X_pca, for the same
# reason 02 rebuilt it: the atlas-level PCA is dominated by whole-kidney
# cell-type differences, not substructure within 715 Uro cells.
#
# Cluster labels for Slingshot's start.clus: FindClusters(graph.name="wsnn")
# (FindMultiModalNeighbors' default SNN graph name -- not "wnn") produces
# NEW cluster IDs (not the RNA-only seurat_clusters from 02), so
# the E16.5-enriched root cluster is picked programmatically here (by %E16.5
# composition) rather than hardcoding an ID the way 03 hardcodes start.clus
# = "4" from a one-time visual inspection of Figure 4.
#
# Input:  UrotheliumDevelopmentScripts/output/UrotheliumOnly_RNA_ATAC.rds
#         (written by 08_UrotheliumDevelopment_ATAC_Build.R)
#         UrotheliumDevelopmentScripts/output/UrotheliumOnly_pseudotime.rds
#         (written by 03_UrotheliumDevelopment_Pseudotime.R -- RNA-only
#         Slingshot pseudotime, for the Panel B comparison)
# Output: FigA_UrotheliumDevelopment_WNN_UMAP_ByStage.pdf
#         FigB_UrotheliumDevelopment_WNN_vs_RNA_Pseudotime_Scatter.pdf
#         UrotheliumOnly_RNA_ATAC_WNN.rds
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(slingshot)
  library(ggplot2)
  library(patchwork)
  library(viridisLite)
  library(dplyr)
})

SCRIPT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/UrotheliumDevelopmentScripts"
OUT_DIR    <- file.path(SCRIPT_DIR, "output")

STAGE_ORDER  <- c("E16.5", "P0", "W3", "W12", "W52", "W92")
STAGE_COLORS <- setNames(viridisLite::viridis(length(STAGE_ORDER)), STAGE_ORDER)

# ── Load ─────────────────────────────────────────────────────────────────────
message("==> Loading combined RNA+ATAC Uro object ...")
in_rds <- file.path(OUT_DIR, "UrotheliumOnly_RNA_ATAC.rds")
if (!file.exists(in_rds)) {
  stop("Missing ", in_rds, " -- run 08_UrotheliumDevelopment_ATAC_Build.R first.")
}
uro <- readRDS(in_rds)
uro$Age <- factor(uro$Age, levels = STAGE_ORDER)
stopifnot("pca" %in% Reductions(uro))
message(sprintf("  %d cells, ATAC: %d peaks, RNA pca dims: %d",
                 ncol(uro), nrow(uro[["ATAC"]]), ncol(Embeddings(uro, "pca"))))

################################################################################
# ATAC LSI
################################################################################
message("\n==> ATAC dimensionality reduction (TF-IDF -> SVD) ...")
DefaultAssay(uro) <- "ATAC"
uro <- RunTFIDF(uro)
uro <- FindTopFeatures(uro, min.cutoff = "q0")
uro <- RunSVD(uro)

# Component 1 is near-universally a sequencing-depth artifact (standard
# Signac practice) -- confirm via correlation with total ATAC counts before
# dropping it, rather than dropping blindly.
depth_cor <- cor(Embeddings(uro, "lsi")[, 1], uro$atac_total_counts, method = "spearman")
message(sprintf(
  "  LSI component 1 vs. total ATAC counts: Spearman rho = %.3f (dropped from WNN regardless, per standard practice)",
  depth_cor))
lsi_dims <- 2:min(30, ncol(Embeddings(uro, "lsi")))

################################################################################
# WNN integration
################################################################################
message("\n==> FindMultiModalNeighbors (RNA pca + ATAC lsi) ...")
uro <- FindMultiModalNeighbors(
  uro,
  reduction.list = list("pca", "lsi"),
  dims.list      = list(1:ncol(Embeddings(uro, "pca")), lsi_dims)
)

# 2D embedding for Panel A / plotting.
uro <- RunUMAP(uro, nn.name = "weighted.nn", reduction.name = "wnn.umap",
               reduction.key = "wnnUMAP_")
# Separate higher-dim embedding purely for Slingshot fitting (see header).
uro <- RunUMAP(uro, nn.name = "weighted.nn", reduction.name = "wnn.umap.fit",
               reduction.key = "wnnUMAPfit_", n.components = 10)
uro <- FindClusters(uro, graph.name = "wsnn", resolution = 0.5, verbose = FALSE)
uro$wnn_clusters <- uro$seurat_clusters

message("  WNN cluster x stage crosstab:")
print(table(uro$wnn_clusters, uro$Age))

################################################################################
# Figure A: WNN UMAP by stage
################################################################################
message("\n==> Building Figure A (WNN UMAP by stage) ...")

figA <- DimPlot(uro, reduction = "wnn.umap", group.by = "Age",
                 cols = STAGE_COLORS, pt.size = 1.3) +
  ggtitle("Panel A. Joint WNN UMAP, colored by developmental stage") +
  labs(color = "Stage") +
  coord_fixed() +
  theme(aspect.ratio = 1, plot.title = element_text(face = "bold", size = 13))

ggsave(file.path(OUT_DIR, "FigA_UrotheliumDevelopment_WNN_UMAP_ByStage.pdf"),
       figA, width = 7, height = 6)
message("  Saved: FigA_UrotheliumDevelopment_WNN_UMAP_ByStage.pdf")

################################################################################
# Panel B: WNN-based Slingshot pseudotime vs. RNA-only Slingshot pseudotime
################################################################################
message("\n==> Fitting Slingshot on the 10-dim WNN embedding ...")

# Root: whichever WNN cluster has the highest %E16.5 composition (mirrors
# 03's "most purely-fetal subcluster" root logic, but computed programmatically
# since these are new WNN cluster IDs, not 02/03's RNA-only subcluster IDs).
stage_frac <- table(uro$wnn_clusters, uro$Age == "E16.5")
pct_e165 <- stage_frac[, "TRUE"] / rowSums(stage_frac)
root_clus <- names(which.max(pct_e165))
message(sprintf("  Root cluster: %s (%.0f%% E16.5, n=%d)", root_clus,
                 100 * pct_e165[root_clus], sum(uro$wnn_clusters == root_clus)))

wnn_fit_emb <- Embeddings(uro, "wnn.umap.fit")
sds_wnn <- slingshot(wnn_fit_emb, clusterLabels = uro$wnn_clusters, start.clus = root_clus)
pt_mat_wnn <- slingPseudotime(sds_wnn)
uro$WNN_Pseudotime <- rowMeans(pt_mat_wnn, na.rm = TRUE)

message("  WNN pseudotime range: ", paste(round(range(uro$WNN_Pseudotime, na.rm = TRUE), 2), collapse = " - "))
message("  Mean WNN pseudotime by stage:")
print(uro@meta.data %>% group_by(Age) %>% summarise(mean_pt = mean(WNN_Pseudotime, na.rm = TRUE), .groups = "drop"))

message("\n==> Loading RNA-only Slingshot pseudotime for comparison ...")
rna_pt_rds <- file.path(OUT_DIR, "UrotheliumOnly_pseudotime.rds")
if (!file.exists(rna_pt_rds)) {
  stop("Missing ", rna_pt_rds, " -- run 03_UrotheliumDevelopment_Pseudotime.R first.")
}
rna_pt_obj <- readRDS(rna_pt_rds)
stopifnot(identical(colnames(rna_pt_obj), colnames(uro)) || setequal(colnames(rna_pt_obj), colnames(uro)))
uro$RNA_Pseudotime <- rna_pt_obj$Pseudotime[match(colnames(uro), colnames(rna_pt_obj))]

pt_df <- uro@meta.data %>%
  filter(!is.na(WNN_Pseudotime), !is.na(RNA_Pseudotime)) %>%
  select(Age, WNN_Pseudotime, RNA_Pseudotime)
rho_wnn_rna <- cor(pt_df$WNN_Pseudotime, pt_df$RNA_Pseudotime, method = "spearman")
message(sprintf("\n  Spearman correlation(WNN pseudotime, RNA-only Slingshot pseudotime) = %.3f (n=%d cells)",
                 rho_wnn_rna, nrow(pt_df)))

saveRDS(uro, file.path(OUT_DIR, "UrotheliumOnly_RNA_ATAC_WNN.rds"))
message("  Saved: UrotheliumOnly_RNA_ATAC_WNN.rds")

################################################################################
# Figure B: WNN vs. RNA pseudotime scatter
################################################################################
message("\n==> Building Figure B (WNN vs. RNA pseudotime scatter) ...")

figB <- ggplot(pt_df, aes(x = RNA_Pseudotime, y = WNN_Pseudotime, color = Age)) +
  geom_point(size = 1.6, alpha = 0.75) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.6, se = TRUE, inherit.aes = FALSE,
              aes(x = RNA_Pseudotime, y = WNN_Pseudotime)) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4,
           label = sprintf("rho = %.2f", rho_wnn_rna), size = 4) +
  scale_color_manual(values = STAGE_COLORS) +
  labs(x = "RNA-only Slingshot pseudotime (script 03)",
       y = "WNN Slingshot pseudotime (this script)", color = "Stage",
       title = "Panel B. WNN-based vs. RNA-only pseudotime",
       subtitle = "Cell-by-cell comparison, matched Uro cells") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13))

ggsave(file.path(OUT_DIR, "FigB_UrotheliumDevelopment_WNN_vs_RNA_Pseudotime_Scatter.pdf"),
       figB, width = 7, height = 6)
message("  Saved: FigB_UrotheliumDevelopment_WNN_vs_RNA_Pseudotime_Scatter.pdf")

message("\n==> Done.")
