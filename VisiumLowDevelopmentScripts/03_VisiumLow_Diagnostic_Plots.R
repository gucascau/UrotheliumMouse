################################################################################
# 03_VisiumLow_Diagnostic_Plots.R
#
# Diagnostic plots for the Harmony-integrated Visium developmental-kidney
# object (02_VisiumLow_Integrate_Harmony.R): UMAP by cluster/stage/sample,
# plus a per-sample spatial cluster map (24 samples). Split out from the
# integration script so plots can be regenerated without repeating the
# load/merge/Harmony/cluster step.
#
# Spatial panels are built manually with ggplot2 (annotation_custom +
# geom_point) rather than Seurat's SpatialDimPlot: SpatialDimPlot renders a
# blank/gray gradient instead of the tissue image for these particular
# objects (verified the raw @image array itself is valid RGB data via base R
# plot() -- the bug is specifically in Seurat's ggplot rendering path for
# these updated-from-old-SeuratObject VisiumV1 images, not the data). The
# manual approach reads @image + @coordinates directly and bypasses it.
#
# Input:  VisiumLowDevelopmentScripts/output/VisiumLow_harmony_integrated.rds
# Output: VisiumLow_harmony_clusters.pdf
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(viridisLite)
  library(grid)
})

# ── Manual spatial plot helper (bypasses SpatialDimPlot's broken image render) ──
plot_spatial_manual <- function(obj, image_name, color_var, title = "",
                                 color_values = NULL, continuous = FALSE,
                                 legend_title = color_var, pt_size = 1.2) {
  img_s4 <- obj@images[[image_name]]
  coords <- img_s4@coordinates
  coords$color_val <- obj[[color_var, drop = TRUE]][rownames(coords)]
  img_arr <- img_s4@image
  img_dim <- dim(img_arr)
  g <- rasterGrob(img_arr, width = unit(1, "npc"), height = unit(1, "npc"))

  p <- ggplot(coords, aes(x = imagecol, y = imagerow)) +
    annotation_custom(g, xmin = 0, xmax = img_dim[2], ymin = 0, ymax = img_dim[1]) +
    geom_point(aes(color = color_val), size = pt_size) +
    scale_y_reverse() +
    coord_fixed(xlim = c(0, img_dim[2]), ylim = c(img_dim[1], 0)) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA)) +
    ggtitle(title)

  if (continuous) {
    p <- p + scale_color_viridis_c(option = "plasma", name = legend_title)
  } else if (!is.null(color_values)) {
    p <- p + scale_color_manual(values = color_values, name = legend_title)
  } else {
    p <- p + labs(color = legend_title)
  }
  p
}

SCRIPT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/VisiumLowDevelopmentScripts"
OUT_DIR    <- file.path(SCRIPT_DIR, "output")

STAGE_ORDER  <- c("E16.5", "P0", "W3", "W12", "W52", "W92")
STAGE_COLORS <- setNames(viridisLite::viridis(length(STAGE_ORDER)), STAGE_ORDER)

message("==> Loading integrated object ...")
in_rds <- file.path(OUT_DIR, "VisiumLow_harmony_integrated.rds")
if (!file.exists(in_rds)) {
  stop("Missing ", in_rds, " -- run 02_VisiumLow_Integrate_Harmony.R first.")
}
object <- readRDS(in_rds)
object$Age <- factor(object$Age, levels = STAGE_ORDER)
message(sprintf("  %d spots, %d samples, %d clusters", ncol(object),
                length(unique(object$sample_id)), length(unique(object$seurat_cluster.harmony))))

out_pdf <- file.path(OUT_DIR, "VisiumLow_harmony_clusters.pdf")
pdf(out_pdf, width = 14, height = 6)

message("==> UMAP overview plots ...")
# shuffle = TRUE: cells were merged in Age order (01_VisiumLow_BuildMetadata.R
# sorts the sample table by Age), so without shuffling, ggplot draws points in
# that order and the last-loaded stage (W92) visually overplots everything
# else regardless of its actual share of spots -- shuffle removes that
# z-order artifact.
p_cluster <- DimPlot(object, reduction = "umap.harmony", shuffle = TRUE, seed = 1,
                      group.by = "seurat_cluster.harmony", label = TRUE, repel = TRUE) +
  ggtitle("Harmony clusters")
p_stage <- DimPlot(object, reduction = "umap.harmony", shuffle = TRUE, seed = 1,
                    group.by = "Age", cols = STAGE_COLORS) +
  ggtitle("Developmental stage")
print(p_cluster | p_stage)

p_sample <- DimPlot(object, reduction = "umap.harmony", group.by = "sample_id",
                     shuffle = TRUE, seed = 1) +
  ggtitle("Sample identity (batch-correction check)") +
  theme(legend.text = element_text(size = 6)) +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))
print(p_sample)

message("==> Per-sample spatial cluster maps (24 samples) ...")
for (samp in unique(object$sample_id)) {
  message("  ", samp)
  sub_obj <- suppressWarnings(subset(object, cells = colnames(object)[object$sample_id == samp]))
  img_name <- Images(sub_obj)[1]
  tryCatch({
    print(
      plot_spatial_manual(sub_obj, img_name, "seurat_cluster.harmony",
                           title = samp, legend_title = "Cluster", pt_size = 1.2)
    )
  }, error = function(e) {
    message("    Spatial plot failed for ", samp, ": ", e$message)
  })
}

dev.off()
message("\n  Saved: ", out_pdf)
message("\n==> Done.")
