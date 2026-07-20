################################################################################
# 03_UrotheliumDevelopment_Pseudotime.R
#
# Trajectory / pseudotime analysis of the reclustered Uro-only population
# (02_UrotheliumDevelopment_Subcluster_UMAP.R), using slingshot: a
# cluster-based minimum spanning tree fit on the Uro-only PCA, principal
# curves per lineage, and per-cell pseudotime along those curves. Curves are
# then re-embedded into the Uro-only UMAP purely for visualization.
#
# Figure 5. Developmental trajectory
#   UMAP colored by (A) stage and (B) pseudotime, with the fitted slingshot
#   curve(s) overlaid, arrow indicating the direction of increasing
#   pseudotime (rooted at the most fetal subcluster). (C) Pseudotime
#   distribution by stage, with the Spearman correlation annotated -- makes
#   directly visible what that correlation number means: a clean fetal vs.
#   postnatal separation with a flat, overlapping postnatal tail.
#
# Figure 6. Smoothed marker expression along pseudotime
#   6 representative genes (2 per marker set from Figure 2: upper-tract
#   developmental, immature/proliferative, differentiation/barrier),
#   scatter of expression vs. pseudotime (points colored by stage) with a
#   GAM-smoothed trend per gene.
#
# Root: start.clus = "4", the most purely-fetal subcluster from Figure 4
# (38/42 cells are E16.5) -- this fixes the trajectory's direction to run
# fetal -> mature rather than leaving slingshot to pick an arbitrary root.
#
# Input:  UrotheliumDevelopmentScripts/output/UrotheliumOnly_reclustered.rds
#         (written by 02_UrotheliumDevelopment_Subcluster_UMAP.R)
# Output: Fig5_UrotheliumDevelopment_Trajectory_UMAP.pdf
#         Fig6_UrotheliumDevelopment_SmoothedExpression_Pseudotime.pdf
#         UrotheliumOnly_pseudotime.rds
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(slingshot)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(viridisLite)
})

SCRIPT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/UrotheliumDevelopmentScripts"
OUT_DIR    <- file.path(SCRIPT_DIR, "output")

STAGE_ORDER  <- c("E16.5", "P0", "W3", "W12", "W52", "W92")
STAGE_COLORS <- setNames(viridisLite::viridis(length(STAGE_ORDER)), STAGE_ORDER)

# ── Load the reclustered Uro-only object ────────────────────────────────────
message("==> Loading reclustered Uro-only object ...")
recl_rds <- file.path(OUT_DIR, "UrotheliumOnly_reclustered.rds")
if (!file.exists(recl_rds)) {
  stop("Missing ", recl_rds, " -- run 02_UrotheliumDevelopment_Subcluster_UMAP.R first.")
}
uro <- readRDS(recl_rds)
uro$Age <- factor(uro$Age, levels = STAGE_ORDER)
message(sprintf("  %d cells, %d subclusters", ncol(uro), length(levels(uro$seurat_clusters))))

################################################################################
# Slingshot trajectory
################################################################################
message("\n==> Fitting slingshot trajectory (PCA space, root = subcluster 4) ...")

pca_emb  <- Embeddings(uro, "pca")
umap_emb <- Embeddings(uro, "umap")

sds <- slingshot(pca_emb, clusterLabels = uro$seurat_clusters, start.clus = "4")

lineages <- slingLineages(sds)
message(sprintf("  Number of lineages detected: %d", length(lineages)))
for (i in seq_along(lineages)) {
  message(sprintf("    Lineage %d: %s", i, paste(lineages[[i]], collapse = " -> ")))
}
message("  Cluster-level minimum spanning tree:")
print(slingMST(sds))
if (length(lineages) == 1) {
  message("  -> Single lineage: the subclusters form one linear path with no",
          " branch point (fetal subclusters feed into one shared mature state,",
          " not multiple diverging fates).")
} else {
  message("  -> Multiple lineages: the trajectory branches into distinct fates",
          " downstream of a shared progenitor-like state.")
}

# Re-embed the fitted curve(s) into the Uro-only UMAP purely for plotting.
sds_umap <- embedCurves(sds, umap_emb)
curve_df <- slingCurves(sds_umap, as.df = TRUE) %>% arrange(Lineage, Order)

# Per-cell pseudotime: average across lineages a cell belongs to (only one
# lineage here, so this is just that lineage's pseudotime; written generally
# in case a future re-run with a different resolution produces branches).
pt_mat <- slingPseudotime(sds)
uro$Pseudotime <- rowMeans(pt_mat, na.rm = TRUE)

message("\n  Pseudotime range: ", paste(round(range(uro$Pseudotime), 2), collapse = " - "))
message("  Mean pseudotime by stage (sanity check -- should increase E16.5 -> W92):")
print(uro@meta.data %>% group_by(Age) %>% summarise(mean_pt = mean(Pseudotime), .groups = "drop"))
stage_rank <- as.integer(uro$Age)
rho_stage <- cor(uro$Pseudotime, stage_rank, method = "spearman")
message(sprintf("  Spearman correlation(pseudotime, chronological stage) = %.3f", rho_stage))

saveRDS(uro, file.path(OUT_DIR, "UrotheliumOnly_pseudotime.rds"))
message("  Saved: UrotheliumOnly_pseudotime.rds")

################################################################################
# Figure 5: trajectory UMAP
################################################################################
message("\n==> Building Figure 5 (trajectory UMAP) ...")

arrow_spec <- arrow(length = unit(0.12, "inches"), type = "closed")

p5a <- DimPlot(uro, reduction = "umap", group.by = "Age", cols = STAGE_COLORS,
               pt.size = 1.3) +
  geom_path(data = curve_df, aes(x = umap_1, y = umap_2, group = Lineage),
            color = "black", linewidth = 0.8, arrow = arrow_spec) +
  ggtitle("Stage + fitted trajectory") +
  labs(color = "Stage") +
  coord_fixed() +
  theme(aspect.ratio = 1, plot.title = element_text(face = "bold"))

umap_df <- data.frame(umap_emb, Pseudotime = uro$Pseudotime)

p5b <- ggplot(umap_df, aes(x = umap_1, y = umap_2)) +
  geom_point(aes(color = Pseudotime), size = 1.3) +
  geom_path(data = curve_df, aes(x = umap_1, y = umap_2, group = Lineage),
            color = "black", linewidth = 0.8, arrow = arrow_spec) +
  scale_color_viridis_c(option = "plasma") +
  labs(x = "umap_1", y = "umap_2") +
  ggtitle("Pseudotime + fitted trajectory") +
  coord_fixed() +
  theme_classic() +
  theme(aspect.ratio = 1, plot.title = element_text(face = "bold"),
        legend.title = element_text(size = 9))

p5c <- ggplot(uro@meta.data, aes(x = Age, y = Pseudotime, fill = Age)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
  scale_fill_manual(values = STAGE_COLORS, guide = "none") +
  labs(x = "Developmental stage", y = "Pseudotime",
       subtitle = sprintf("Spearman rho = %.2f (vs. chronological stage)", rho_stage)) +
  ggtitle("Pseudotime by stage") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

fig5 <- (p5a | p5b) / p5c +
  plot_annotation(title = "Figure 5. Urothelium developmental trajectory",
                   theme = theme(plot.title = element_text(face = "bold", size = 14)))

ggsave(file.path(OUT_DIR, "Fig5_UrotheliumDevelopment_Trajectory_UMAP.pdf"),
       fig5, width = 12, height = 10)
message("  Saved: Fig5_UrotheliumDevelopment_Trajectory_UMAP.pdf")

################################################################################
# Figure 6: smoothed expression along pseudotime
################################################################################
message("\n==> Building Figure 6 (smoothed expression vs. pseudotime) ...")

# 2 genes from each Figure-2 marker set, spanning the fetal -> mature program.
PSEUDOTIME_GENES <- c("Pkhd1", "Bicc1",   # upper-tract developmental
                       "Trp63", "Mki67",  # immature / proliferative
                       "Upk1b", "Krt20")  # differentiation / barrier

fm  <- uro[["originalexp"]]@meta.features
idx <- match(PSEUDOTIME_GENES, fm$gene_symbols)
if (any(is.na(idx))) {
  stop("Marker symbols not found in gene_symbols: ",
       paste(PSEUDOTIME_GENES[is.na(idx)], collapse = ", "))
}
ens_ids <- rownames(fm)[idx]

expr_mat <- GetAssayData(uro, assay = "originalexp", layer = "data")[ens_ids, , drop = FALSE]
rownames(expr_mat) <- PSEUDOTIME_GENES

expr_long <- as.data.frame(t(as.matrix(expr_mat))) %>%
  tibble::rownames_to_column("cell") %>%
  pivot_longer(-cell, names_to = "Gene", values_to = "Expression") %>%
  mutate(Gene = factor(Gene, levels = PSEUDOTIME_GENES)) %>%
  left_join(
    data.frame(cell = colnames(uro), Pseudotime = uro$Pseudotime, Age = uro$Age),
    by = "cell"
  )

fig6 <- ggplot(expr_long, aes(x = Pseudotime, y = Expression)) +
  geom_point(aes(color = Age), size = 0.9, alpha = 0.5) +
  # k = 4: with n = 715 and heavy dropout/zero-inflation, the GAM default
  # (k = 10) fits small wiggles that are dropout noise, not biology; a lower
  # basis dimension keeps only the broad rise/fall/plateau trend.
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 4),
              color = "black", linewidth = 0.8, se = TRUE) +
  scale_color_manual(values = STAGE_COLORS) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 3) +
  labs(x = "Pseudotime", y = "Log-normalized expression", color = "Stage",
       title = "Figure 6. Smoothed marker expression along pseudotime") +
  theme_classic(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 14)) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

ggsave(file.path(OUT_DIR, "Fig6_UrotheliumDevelopment_SmoothedExpression_Pseudotime.pdf"),
       fig6, width = 10, height = 6)
message("  Saved: Fig6_UrotheliumDevelopment_SmoothedExpression_Pseudotime.pdf")

message("\n==> Done.")
