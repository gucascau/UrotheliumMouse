################################################################################
# 02c_UMAP_plots_AllscRNA.R
#
# Input : output/AllUrothelium_gated_UMAP_nosnRNA.rds
#         (scRNA-seq only urothelium cells, Harmony-integrated)
#
# Outputs (output/AllUrothelium_nosnRNA_plots/):
#   1. FeaturePlot  — urothelium markers on UMAP
#   2. DotPlot      — markers × each metadata group
#   3. DimPlot      — UMAP coloured / split by each metadata group
#
# Pipeline:
#   STEP 8 : FindNeighbors + FindClusters + RunUMAP on Harmony embedding
#   STEP 9 : Marker & group visualisations
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})
options(future.globals.maxSize = 8 * 1024^3)

# ── Parameters ────────────────────────────────────────────────────────────────
HARMONY_DIMS <- 1:50
RESOLUTION   <- 0.5

log_mem <- function(label) {
  g  <- gc(verbose = FALSE)
  mb <- sum(g[, 2])
  message(sprintf("  [mem] %s: %.1f GB in use", label, mb / 1024))
}

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/FinalUrotheliumCells"
SCR_DIR  <- file.path(BASE_DIR, "UrotheliumScripts")
OUT_DIR  <- file.path(SCR_DIR, "output")
#PLOT_DIR <- file.path(OUT_DIR, "AllUrothelium_nosnRNA_plots")
PLOT_DIR <- file.path(OUT_DIR, "AllUrothelium_withsnRNA_plots")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

args     <- commandArgs(trailingOnly = TRUE)
RDS_PATH <- if (length(args) >= 1) args[[1]] else
  file.path(OUT_DIR, "AllUrothelium_gated_UMAP.rds")

# ── Marker genes ──────────────────────────────────────────────────────────────
URO_MARKERS <- c(
  "Krt8",  "Krt18", "Krt19",
  "Upk1a", "Upk1b", "Upk2", "Upk3a", "Upk3b",
  "Krt20", "Krt5",  "Krt14", "Trp63",
  "Foxa1", "Gata3", "Pparg"
)

DOTPLOT_GROUPS <- c("seurat_clusters",  "FinalConditionL1", "FinalConditionL2", "Finalgsm_id", "FinalSampleId", "Finaltechnology", "Finalpaper", "Finalsource_GEO", "Finalgsm_id", "Finaltissue")


# ── Load ──────────────────────────────────────────────────────────────────────
message("Loading: ", RDS_PATH)
if (!file.exists(RDS_PATH))
  stop("Input RDS not found: ", RDS_PATH)

so <- readRDS(RDS_PATH)
message(sprintf("  Loaded: %d cells x %d genes", ncol(so), nrow(so)))
log_mem("after load")

DimPlot(so, label = TRUE, repel = TRUE)

# Restore data layer if needed
rna_layers <- Layers(so[["RNA"]])
if (!"data" %in% rna_layers) {
  message("  data layer missing -- joining layers ...")
  so <- JoinLayers(so)
}

# Verify markers
present <- intersect(URO_MARKERS, rownames(so))
missing <- setdiff(URO_MARKERS, rownames(so))
if (length(missing) > 0)
  message(sprintf("  Markers not found (skipped): %s",
                  paste(missing, collapse = ", ")))
if (length(present) == 0)
  stop("None of the requested marker genes found in the object.")
message(sprintf("  Plotting %d / %d markers: %s",
                length(present), length(URO_MARKERS),
                paste(present, collapse = ", ")))




################################################################################
# STEP 9: FeaturePlot -- all markers on UMAP
################################################################################

message("Generating FeaturePlots ...")
n_col <- 5
n_row <- ceiling(length(present) / n_col)

fp <- FeaturePlot(
  so,
  features  = present,
  reduction = "umap_harmony",
  ncol      = n_col,
  raster    = TRUE,
  order     = TRUE,
  cols      = c("lightgrey", "#d73027")
)
pdf(file.path(PLOT_DIR, "AllUrothelium_FeaturePlot_markers.pdf"),
    width = n_col * 4, height = n_row * 4)
print(fp)
dev.off()
message("  Saved: AllUrothelium_FeaturePlot_markers.pdf")

fp_split <- FeaturePlot(
  so,
  features  = present,
  reduction = "umap_harmony",
  split.by  = "Categories",
  ncol      = n_col,
  raster    = TRUE,
  order     = TRUE,
  cols      = c("lightgrey", "#d73027")
)
pdf(file.path(PLOT_DIR, "AllUrothelium_FeaturePlot_markers_splitCategories.pdf"),
    width = n_col * 8, height = n_row * 4)
print(fp_split)
dev.off()
message("  Saved: AllUrothelium_FeaturePlot_markers_splitCategories.pdf")

# Individual FeaturePlot PDFs (one per gene)
message("Generating individual gene FeaturePlots ...")
for (gene in present) {
  p <- FeaturePlot(
    so,
    features  = gene,
    reduction = "umap_harmony",
    raster    = TRUE,
    order     = TRUE,
    cols      = c("lightgrey", "#d73027")
  ) + ggtitle(gene)
  ggsave(
    file.path(PLOT_DIR, sprintf("FeaturePlot_%s.pdf", gene)),
    plot = p, width = 6, height = 5
  )
}
message(sprintf("  Saved %d individual FeaturePlot PDFs", length(present)))


################################################################################
# STEP 9: DotPlot -- markers x each metadata grouping
################################################################################

message("Generating DotPlots ...")
for (grp in DOTPLOT_GROUPS) {
  if (!grp %in% colnames(so@meta.data)) next

  grp_vals <- as.character(so@meta.data[[grp]])
  grp_vals[is.na(grp_vals) | grp_vals == ""] <- NA
  valid_levels <- names(which(table(grp_vals) >= 3))
  if (length(valid_levels) == 0) {
    message(sprintf("  Skipping DotPlot for %s: no groups with >= 3 cells", grp))
    next
  }
  so_sub <- so[, !is.na(grp_vals) & grp_vals %in% valid_levels]
  so_sub@meta.data[[grp]] <- droplevels(
    factor(so_sub@meta.data[[grp]], levels = valid_levels)
  )
  n_dropped <- ncol(so) - ncol(so_sub)
  if (n_dropped > 0)
    message(sprintf("  DotPlot %s: dropped %d cells in sparse/NA groups",
                    grp, n_dropped))

  n_levels <- length(valid_levels)
  height   <- max(4, 2 + n_levels * 0.35)
  width    <- max(8, 2 + length(present) * 0.6)

  p <- DotPlot(
    so_sub,
    features  = present,
    group.by  = grp,
    cols      = c("lightgrey", "#08519c"),
    dot.scale = 6,
    scale     = TRUE
  ) +
    scale_color_gradientn(
      colors   = c("lightgrey", "#08519c"),
      na.value = "lightgrey"     # NaN from all-zero groups → lightgrey, not grey50
    ) +
    RotatedAxis() +
    labs(title = sprintf("Urothelium markers -- grouped by %s", grp)) +
    theme(axis.text.x = element_text(size = 9),
          axis.text.y = element_text(size = 9))

  fname <- sprintf("AllUrothelium_DotPlot_by_%s.pdf", grp)
  ggsave(file.path(PLOT_DIR, fname), plot = p, width = width, height = height)
  message(sprintf("  Saved: %s", fname))
}

# Combined: DimPlot clusters + DotPlot by cluster
message("Generating combined overview plot ...")
p_dim <- DimPlot(
  so, group.by = "seurat_clusters", reduction = "umap_harmony",
  label = TRUE, repel = TRUE, raster = TRUE
) + ggtitle("Clusters") + NoLegend()

p_dot <- DotPlot(
  so, features = present, group.by = "seurat_clusters",
  cols = c("lightgrey", "#08519c"), dot.scale = 5, scale = TRUE
) +
  scale_color_gradientn(colors = c("lightgrey", "#08519c"), na.value = "lightgrey") +
  RotatedAxis() +
  labs(title = "Markers by cluster") +
  theme(axis.text.x = element_text(size = 8))

pdf(file.path(PLOT_DIR, "AllUrothelium_overview_clusters_dotplot.pdf"),
    width = 22, height = 10)
print(p_dim | p_dot)
dev.off()
message("  Saved: AllUrothelium_overview_clusters_dotplot.pdf")

# DotPlot panels split by tissue
if ("tissue" %in% colnames(so@meta.data)) {
  tissues <- sort(unique(so@meta.data$tissue))
  dot_by_tissue <- lapply(tissues, function(tis) {
    so_sub <- so[, so$tissue == tis]
    if (ncol(so_sub) < 10) return(NULL)
    DotPlot(so_sub, features = present, group.by = "seurat_clusters",
            cols = c("lightgrey", "#08519c"), dot.scale = 5, scale = TRUE) +
      scale_color_gradientn(colors = c("lightgrey", "#08519c"), na.value = "lightgrey") +
      RotatedAxis() +
      ggtitle(sprintf("Tissue: %s", tis)) +
      theme(axis.text.x = element_text(size = 7))
  })
  dot_by_tissue <- Filter(Negate(is.null), dot_by_tissue)
  if (length(dot_by_tissue) > 0) {
    pdf(file.path(PLOT_DIR, "AllUrothelium_DotPlot_byTissue_panels.pdf"),
        width = 20, height = 8 * length(dot_by_tissue))
    print(wrap_plots(dot_by_tissue, ncol = 1))
    dev.off()
    message("  Saved: AllUrothelium_DotPlot_byTissue_panels.pdf")
  }
}


################################################################################
# STEP 9: DimPlot -- UMAP coloured / split by each metadata group
################################################################################

message("Generating DimPlots ...")

so@meta.data %>% head()

# analyses for FinalConditionL1, FinalConditionL2, Finalgsm_id, FinalSampleId, Finaltechnology, Finalpaper, Finalsource_GEO, Finalgsm_id, Finaltissue
DOTPLOT_GROUPS 
for (grp in DOTPLOT_GROUPS) {
  # we mannually adjusted the condition figure size to better fit the legend, so skip it here
  # we manually ajusted the sample_id figure size to better fit the legend, so skip it here
  # grp = c("FinalSampleId")
  if (!grp %in% colnames(so@meta.data)) next
  n_levels <- length(unique(na.omit(so@meta.data[[grp]])))
  dim_size <- max(7, 3 + n_levels * 0.1)

  # 
  p <- DimPlot(
    so, group.by = grp, reduction = "umap_harmony",
    label = FALSE, repel = TRUE, raster = TRUE
  ) +
    labs(title = sprintf("UMAP -- grouped by %s", grp)) +
    theme(axis.text.x = element_text(size = 9),
          axis.text.y = element_text(size = 9))

  fname <- sprintf("AllUrothelium_DimPlot_by_%s.pdf", grp)
  ggsave(file.path(PLOT_DIR, fname), plot = p,
         width = dim_size, height = dim_size)
  
  # ggsave(file.path(PLOT_DIR, fname), plot = p, width = 20, height = 7)
  message(sprintf("  Saved: %s", fname))
}

for (grp in DOTPLOT_GROUPS) {
  if (!grp %in% colnames(so@meta.data)) next
  n_levels <- length(unique(na.omit(so@meta.data[[grp]])))
  width    <- max(15, 3 + n_levels * 0.6)
  height   <- max(5, 3 + 0.35)

  p <- DimPlot(
    so, group.by = "seurat_clusters", split.by = grp,
    label = FALSE, repel = TRUE, raster = TRUE
  ) +
    NoLegend() +
    labs(title = sprintf("UMAP -- split by %s", grp)) +
    theme(axis.text.x = element_text(size = 9),
          axis.text.y = element_text(size = 9))

  fname <- sprintf("AllUrothelium_DimPlot_Splitby_%s.pdf", grp)
  ggsave(file.path(PLOT_DIR, fname), plot = p, width = width, height = height)
  message(sprintf("  Saved: %s", fname))
}

message("\nAll plots saved to: ", PLOT_DIR)
