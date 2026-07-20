################################################################################
# 02_UrotheliumDevelopment_Subcluster_UMAP.R
#
# Recluster the urothelial (Uro) population from the Chen2025 developmental
# atlas on its own variable genes/PCA. The atlas-level X_pca/X_umap
# (01_UrotheliumDevelopment_Figures.R) were computed across all 203,139
# whole-kidney cells and are dominated by major cell-type differences -- they
# are not appropriate for resolving substructure within the small (n = 715)
# Uro population, so PCA/UMAP are recomputed here from Uro-only HVGs.
#
# Only a "data" (log-normalized) layer exists in the source object -- no raw
# counts -- so there is nothing to re-normalize; we go straight from HVGs to
# scaling/PCA on the existing log-normalized values.
#
# Caveat: 715 cells are spread across 32 samples (~22 cells/sample median),
# too few per sample for Harmony-style batch correction to be reliable, so no
# batch integration is attempted here. Sample identity is not regressed out;
# if subclusters turn out to track sample_id rather than biology, that is a
# batch-effect flag worth checking (see the printed subcluster x sample
# crosstab).
#
# Input:  UrotheliumDevelopmentScripts/output/UrotheliumOnly_raw.rds
#         (written by 01_UrotheliumDevelopment_Figures.R)
# Output: Fig4_UrotheliumDevelopment_Subcluster_UMAP.pdf
#         UrotheliumOnly_reclustered.rds
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(viridisLite)
  library(RColorBrewer)
})

SCRIPT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/UrotheliumDevelopmentScripts"
OUT_DIR    <- file.path(SCRIPT_DIR, "output")

STAGE_ORDER  <- c("E16.5", "P0", "W3", "W12", "W52", "W92")
STAGE_COLORS <- setNames(viridisLite::viridis(length(STAGE_ORDER)), STAGE_ORDER)

# ── Load the Uro-only subset ────────────────────────────────────────────────
message("==> Loading Uro-only subset ...")
uro_rds <- file.path(OUT_DIR, "UrotheliumOnly_raw.rds")
if (!file.exists(uro_rds)) {
  stop("Missing ", uro_rds, " -- run 01_UrotheliumDevelopment_Figures.R first.")
}
uro <- readRDS(uro_rds)
uro$Age <- factor(uro$Age, levels = STAGE_ORDER)
DefaultAssay(uro) <- "originalexp"
message(sprintf("  %d Uro cells across %d samples", ncol(uro), length(unique(uro$sample_id))))
print(table(uro$Age))

# ── Recluster on Uro-only variable genes ────────────────────────────────────
set.seed(1)
N_PCS <- 15  # conservative given n = 715 cells

uro <- FindVariableFeatures(uro, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
uro <- ScaleData(uro, verbose = FALSE)
uro <- RunPCA(uro, npcs = N_PCS, verbose = FALSE)
uro <- FindNeighbors(uro, dims = 1:N_PCS, verbose = FALSE)
uro <- FindClusters(uro, resolution = 0.5, verbose = FALSE)
uro <- RunUMAP(uro, dims = 1:N_PCS, verbose = FALSE)

message("  Subcluster sizes:")
print(table(uro$seurat_clusters))
message("  Subcluster x stage crosstab:")
print(table(uro$seurat_clusters, uro$Age))
message("  Subcluster x sample crosstab (batch-effect check):")
print(table(uro$seurat_clusters, uro$sample_id))

saveRDS(uro, file.path(OUT_DIR, "UrotheliumOnly_reclustered.rds"))
message("  Saved: UrotheliumOnly_reclustered.rds")

################################################################################
# Figure 4: reclustered urothelium UMAP
################################################################################
message("\n==> Building Figure 4 (reclustered UMAP) ...")

cluster_levels  <- levels(uro$seurat_clusters)
cluster_colors  <- setNames(
  RColorBrewer::brewer.pal(max(length(cluster_levels), 3), "Set2")[seq_along(cluster_levels)],
  cluster_levels
)

p4a <- DimPlot(uro, reduction = "umap", group.by = "Age",
               cols = STAGE_COLORS, pt.size = 1.5) +
  ggtitle("Developmental stage") +
  labs(color = "Stage") +
  coord_fixed() +
  theme(aspect.ratio = 1, plot.title = element_text(face = "bold"))

p4b <- DimPlot(uro, reduction = "umap", group.by = "seurat_clusters",
               cols = cluster_colors, pt.size = 1.5) +
  ggtitle("Urothelial subclusters") +
  labs(color = "Subcluster") +
  coord_fixed() +
  theme(aspect.ratio = 1, plot.title = element_text(face = "bold"))

fig4 <- (p4a | p4b) +
  plot_annotation(title = "Figure 4. Urothelium subclustering",
                   theme = theme(plot.title = element_text(face = "bold", size = 14)))

ggsave(file.path(OUT_DIR, "Fig4_UrotheliumDevelopment_Subcluster_UMAP.pdf"),
       fig4, width = 11, height = 5)
message("  Saved: Fig4_UrotheliumDevelopment_Subcluster_UMAP.pdf")

message("\n==> Done.")
