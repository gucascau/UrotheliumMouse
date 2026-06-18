################################################################################
# 05_downstream_AllUrothelium_gated.R
#
# Input : output/AllUrothelium_markers_gated.rds
#
# Steps:
#   1. JoinLayers + FindVariableFeatures (3 000 HVGs)
#   2. ScaleData (regress pct_mt)
#   3. RunPCA (20 PCs)
#   4. RunHarmony (sample_id + technology + source)
#   5. FindNeighbors + FindClusters (res = 0.5)
#   6. RunUMAP  → save RDS
#   7. FeaturePlots (combined, individual, split by Categories)
#   8. DotPlots  × each metadata group
#   9. DimPlots  (grouped + split.by each metadata group)
#
# Output: output/AllUrothelium_gated/
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})
options(future.globals.maxSize = 8 * 1024^3)

log_mem <- function(label) {
  g  <- gc(verbose = FALSE)
  mb <- sum(g[, 2])
  message(sprintf("  [mem] %s: %.1f GB in use", label, mb / 1024))
}

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/FinalUrotheliumCells"
SCR_DIR  <- file.path(BASE_DIR, "UrotheliumScripts")
OUT_DIR  <- file.path(SCR_DIR, "output")
PLOT_DIR <- file.path(OUT_DIR, "AllUrothelium_gated")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

args     <- commandArgs(trailingOnly = TRUE)
RDS_IN   <- if (length(args) >= 1) args[[1]] else
  file.path(OUT_DIR, "AllUrothelium_markers_gated.rds")
RDS_OUT  <- file.path(OUT_DIR, "AllUrothelium_gated_UMAP.rds")

# ── Parameters ────────────────────────────────────────────────────────────────
N_HVG        <- 3000
N_PCS        <- 20
HARMONY_DIMS <- 1:20
RESOLUTION   <- 0.5

# ── Marker genes ──────────────────────────────────────────────────────────────
URO_MARKERS <- c(
  "Krt8",  "Krt18", "Krt19",
  "Upk1a", "Upk1b", "Upk2", "Upk3a", "Upk3b",
  "Krt20", "Krt5",  "Krt14", "Trp63",
  "Foxa1", "Gata3", "Pparg"
)

DOTPLOT_GROUPS <- c("seurat_clusters", "condition", "tissue", "sample_id",
                    "technology", "paper", "Categories")

# ── Load ──────────────────────────────────────────────────────────────────────
message("Loading: ", RDS_IN)
if (!file.exists(RDS_IN)) stop("Input RDS not found: ", RDS_IN)

so <- readRDS(RDS_IN)
message(sprintf("  Loaded: %d cells × %d genes", ncol(so), nrow(so)))
log_mem("after load")

colnames(so@meta.data) 
so@meta.data$source %>% table() %>% sort(decreasing = TRUE)
################################################################################
# STEP 1: JoinLayers + FindVariableFeatures
################################################################################

message("\nJoining layers ...")
so <- JoinLayers(so)

message(sprintf("FindVariableFeatures (%d HVGs) ...", N_HVG))
so <- FindVariableFeatures(so, selection.method = "vst",
                           nfeatures = N_HVG, verbose = FALSE)
message(sprintf("  Top 10 HVGs: %s",
                paste(head(VariableFeatures(so), 10), collapse = ", ")))


################################################################################
# STEP 2: ScaleData
################################################################################

if (!"pct_mt" %in% colnames(so@meta.data)) {
  message("Computing pct_mt ...")
  so[["pct_mt"]] <- PercentageFeatureSet(so, pattern = "^mt-|^MT-")
}
n_na <- sum(is.na(so$pct_mt))
if (n_na > 0) {
  so$pct_mt[is.na(so$pct_mt)] <- median(so$pct_mt, na.rm = TRUE)
  message(sprintf("  Imputed %d pct_mt NAs with median", n_na))
}

message("ScaleData (regress pct_mt, HVGs only) ...")
so <- ScaleData(so,
                features        = VariableFeatures(so),
                vars.to.regress = "pct_mt",
                verbose         = FALSE)
log_mem("after ScaleData")


################################################################################
# STEP 3: RunPCA
################################################################################

message(sprintf("RunPCA (%d PCs) ...", N_PCS))
so <- RunPCA(so, npcs = N_PCS, verbose = FALSE)
log_mem("after PCA")


################################################################################
# STEP 4: RunHarmony
################################################################################

harmony_vars <- "sample_id"
for (v in c("technology", "FinalCategories")) {
  if (v %in% colnames(so@meta.data) &&
      length(unique(na.omit(so@meta.data[[v]]))) > 1) {
    harmony_vars <- c(harmony_vars, v)
  }
}

# Fill NAs in every batch variable — Harmony's design-matrix rbind fails on NAs
for (v in harmony_vars) {
  n_na <- sum(is.na(so@meta.data[[v]]))
  if (n_na > 0) {
    message(sprintf("  Filling %d NAs in '%s' → 'Unknown'", n_na, v))
    so@meta.data[[v]][is.na(so@meta.data[[v]])] <- "Unknown"
  }
}

message(sprintf("RunHarmony (batch = %s) ...", paste(harmony_vars, collapse = " + ")))

so <- RunHarmony(
  so,
  group.by.vars  = harmony_vars,
  reduction      = "pca",
  reduction.save = "harmony",
  verbose        = FALSE
)
# Free scale.data and PCA now that Harmony embedding is computed
so[["RNA"]]$scale.data <- NULL
so[["pca"]] <- NULL
gc()
log_mem("after Harmony")


################################################################################
# STEP 5: FindNeighbors + FindClusters
################################################################################

message("FindNeighbors (annoy, k = 20) ...")
so <- FindNeighbors(
  so,
  reduction    = "harmony",
  dims         = HARMONY_DIMS,
  nn.method    = "annoy",
  k.param      = 20,
  annoy.metric = "euclidean",
  n.trees      = 50,
  verbose      = FALSE
)
log_mem("after FindNeighbors")

message(sprintf("FindClusters (resolution = %.2f) ...", RESOLUTION))
so <- FindClusters(so, resolution = RESOLUTION, verbose = FALSE)
message(sprintf("  Clusters: %d", length(unique(so$seurat_clusters))))
gc()


################################################################################
# STEP 6: RunUMAP + Save
################################################################################

message("RunUMAP on Harmony embedding ...")
so <- RunUMAP(so, reduction = "harmony", dims = HARMONY_DIMS,
              reduction.name = "umap_harmony", verbose = FALSE)
log_mem("after UMAP")

message("Saving: ", RDS_OUT)
saveRDS(so, RDS_OUT)
message("  Saved.")


################################################################################
# STEP 7: FeaturePlots
################################################################################

present <- intersect(URO_MARKERS, rownames(so))
if (length(present) == 0) stop("None of the URO_MARKERS found in the object.")
message(sprintf("\n%d / %d markers present: %s",
                length(present), length(URO_MARKERS),
                paste(present, collapse = ", ")))

n_col <- 5
n_row <- ceiling(length(present) / n_col)

message("Generating FeaturePlots ...")

# Combined — all markers on one page
fp_all <- FeaturePlot(
  so, features = present, reduction = "umap_harmony",
  ncol = n_col, raster = TRUE, order = TRUE,
  cols = c("lightgrey", "#d73027")
)
pdf(file.path(PLOT_DIR, "AllUrothelium_FeaturePlot_markers.pdf"),
    width = n_col * 4, height = n_row * 4)
print(fp_all)
dev.off()
message("  Saved: AllUrothelium_FeaturePlot_markers.pdf")

# Split by Categories
fp_split <- FeaturePlot(
  so, features = present, reduction = "umap_harmony",
  split.by = "Categories", ncol = n_col,
  raster = TRUE, order = TRUE,
  cols = c("lightgrey", "#d73027")
)
pdf(file.path(PLOT_DIR, "AllUrothelium_FeaturePlot_markers_splitCategories.pdf"),
    width = n_col * 8, height = n_row * 4)
print(fp_split)
dev.off()
message("  Saved: AllUrothelium_FeaturePlot_markers_splitCategories.pdf")

# Individual — one PDF per gene
for (gene in present) {
  p <- FeaturePlot(
    so, features = gene, reduction = "umap_harmony",
    raster = TRUE, order = TRUE,
    cols = c("lightgrey", "#d73027")
  ) + ggtitle(gene)
  ggsave(file.path(PLOT_DIR, sprintf("FeaturePlot_%s.pdf", gene)),
         plot = p, width = 6, height = 5)
}
message(sprintf("  Saved %d individual FeaturePlot PDFs", length(present)))


################################################################################
# STEP 8: DotPlots × each metadata group
################################################################################

message("\nGenerating DotPlots ...")

for (grp in DOTPLOT_GROUPS) {
  if (!grp %in% colnames(so@meta.data)) next

  grp_vals    <- as.character(so@meta.data[[grp]])
  grp_vals[is.na(grp_vals) | grp_vals == ""] <- NA
  valid_levels <- names(which(table(grp_vals) >= 3))
  if (length(valid_levels) == 0) {
    message(sprintf("  Skipping DotPlot for %s: no groups with >= 3 cells", grp))
    next
  }
  so_sub <- so[, !is.na(grp_vals) & grp_vals %in% valid_levels]
  so_sub@meta.data[[grp]] <- droplevels(factor(so_sub@meta.data[[grp]],
                                               levels = valid_levels))

  n_levels <- length(valid_levels)
  height   <- max(4, 2 + n_levels * 0.4)
  width    <- max(8, 2 + length(present) * 0.6)

  p <- DotPlot(so_sub, features = present, group.by = grp,
               cols = c("lightgrey", "#08519c"), dot.scale = 5) +
    RotatedAxis() +
    labs(title = sprintf("Markers — grouped by %s", grp)) +
    theme(axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 9))

  ggsave(file.path(PLOT_DIR, sprintf("AllUrothelium_DotPlot_by_%s.pdf", grp)),
         plot = p, width = width, height = height)
  message(sprintf("  Saved: AllUrothelium_DotPlot_by_%s.pdf", grp))
}

# DotPlot per tissue (side-by-side panels)
if ("tissue" %in% colnames(so@meta.data)) {
  tissues        <- sort(unique(so@meta.data$tissue))
  dot_by_tissue  <- Filter(Negate(is.null), lapply(tissues, function(tis) {
    so_sub <- so[, so$tissue == tis]
    if (ncol(so_sub) < 10) return(NULL)
    DotPlot(so_sub, features = present, group.by = "seurat_clusters",
            cols = c("lightgrey", "#08519c"), dot.scale = 5) +
      RotatedAxis() +
      ggtitle(sprintf("Tissue: %s", tis)) +
      theme(axis.text.x = element_text(size = 7))
  }))
  if (length(dot_by_tissue) > 0) {
    pdf(file.path(PLOT_DIR, "AllUrothelium_DotPlot_byTissue_panels.pdf"),
        width = 20, height = 8 * length(dot_by_tissue))
    print(wrap_plots(dot_by_tissue, ncol = 1))
    dev.off()
    message("  Saved: AllUrothelium_DotPlot_byTissue_panels.pdf")
  }
}

# Combined overview: DimPlot clusters + DotPlot by cluster
p_dim <- DimPlot(so, group.by = "seurat_clusters", reduction = "umap_harmony",
                 label = TRUE, repel = TRUE, raster = TRUE) +
  ggtitle("Clusters") + NoLegend()
p_dot <- DotPlot(so, features = present, group.by = "seurat_clusters",
                 cols = c("lightgrey", "#08519c"), dot.scale = 5) +
  RotatedAxis() +
  labs(title = "Markers by cluster") +
  theme(axis.text.x = element_text(size = 8))
pdf(file.path(PLOT_DIR, "AllUrothelium_overview_clusters_dotplot.pdf"),
    width = 22, height = 10)
print(p_dim | p_dot)
dev.off()
message("  Saved: AllUrothelium_overview_clusters_dotplot.pdf")


################################################################################
# STEP 9: DimPlots
################################################################################

message("\nGenerating DimPlots ...")

# Grouped by each metadata column
for (grp in DOTPLOT_GROUPS) {
  if (!grp %in% colnames(so@meta.data)) next

  n_levels <- length(unique(so@meta.data[[grp]]))
  width    <- max(8,  2 + n_levels * 0.35)

  p <- DimPlot(so, group.by = grp, reduction = "umap_harmony",
               label = FALSE, repel = TRUE, raster = TRUE) +
    labs(title = sprintf("UMAP — grouped by %s", grp)) +
    theme(axis.text.x = element_text(size = 9),
          axis.text.y = element_text(size = 9))

  ggsave(file.path(PLOT_DIR, sprintf("AllUrothelium_DimPlot_by_%s.pdf", grp)),
         plot = p, width = width, height = width)
  message(sprintf("  Saved: AllUrothelium_DimPlot_by_%s.pdf", grp))
}

# Clusters split by each metadata column
for (grp in DOTPLOT_GROUPS) {
  if (!grp %in% colnames(so@meta.data)) next

  n_levels <- length(unique(so@meta.data[[grp]]))
  width    <- max(15, 2 + length(present) * 0.6)
  height   <- max(4,  2 + n_levels * 0.35)

  p <- DimPlot(so, group.by = "seurat_clusters", split.by = grp,
               label = FALSE, repel = TRUE, raster = TRUE) +
    NoLegend() +
    labs(title = sprintf("UMAP — clusters split by %s", grp)) +
    theme(axis.text.x = element_text(size = 9),
          axis.text.y = element_text(size = 9))

  ggsave(file.path(PLOT_DIR, sprintf("AllUrothelium_DimPlot_Splitby_%s.pdf", grp)),
         plot = p, width = width, height = height)
  message(sprintf("  Saved: AllUrothelium_DimPlot_Splitby_%s.pdf", grp))
}

message("\nAll done. Output: ", PLOT_DIR)
