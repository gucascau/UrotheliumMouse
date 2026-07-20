################################################################################
# 03_VisiumLowDisease_Diagnostic_Plots.R
#
# Diagnostic plots for the Harmony-integrated disease/reference Visium object
# (02_VisiumLowDisease_Integrate_Harmony.R): UMAP by cluster/study/condition/
# sample/chemistry, plus a per-sample spatial cluster map (10 samples). Split
# out from the integration script so plots can be regenerated without
# repeating the load/merge/Harmony/cluster step.
#
# Manual ggplot2 spatial rendering (annotation_custom + geom_point), same
# approach as VisiumLowDevelopmentScripts/03_VisiumLow_Diagnostic_Plots.R --
# reads @image + @coordinates directly rather than SpatialDimPlot.
#
# Input:  VisiumLowDiseaseScripts/output/VisiumLowDisease_harmony_integrated.rds
# Output: VisiumLowDisease_harmony_clusters.pdf
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(viridisLite)
  library(grid)
})

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

SCRIPT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/VisiumLowDiseaseScripts"
OUT_DIR    <- file.path(SCRIPT_DIR, "output")

STUDY_COLORS <- setNames(
  viridisLite::viridis(3, option = "turbo"),
  c("GSE269063", "GSE269884", "10x_Genomics_Datasets")
)

message("==> Loading integrated object ...")
in_rds <- file.path(OUT_DIR, "VisiumLowDisease_harmony_integrated.rds")
if (!file.exists(in_rds)) {
  stop("Missing ", in_rds, " -- run 02_VisiumLowDisease_Integrate_Harmony.R first.")
}
object <- readRDS(in_rds)
message(sprintf("  %d spots, %d samples, %d clusters", ncol(object),
                length(unique(object$sample_id)), length(unique(object$seurat_cluster.harmony))))

out_pdf <- file.path(OUT_DIR, "VisiumLowDisease_harmony_clusters.pdf")
pdf(out_pdf, width = 14, height = 6)

message("==> UMAP overview plots ...")
p_cluster <- DimPlot(object, reduction = "umap.harmony", shuffle = TRUE, seed = 1,
                      group.by = "seurat_cluster.harmony", label = TRUE, repel = TRUE) +
  ggtitle("Harmony clusters")
p_study <- DimPlot(object, reduction = "umap.harmony", shuffle = TRUE, seed = 1,
                    group.by = "study", cols = STUDY_COLORS) +
  ggtitle("Study of origin")
print(p_cluster | p_study)

p_condition <- DimPlot(object, reduction = "umap.harmony", shuffle = TRUE, seed = 1,
                        group.by = "disease_model") +
  ggtitle("Disease model")
p_chem <- DimPlot(object, reduction = "umap.harmony", shuffle = TRUE, seed = 1,
                   group.by = "chemistry") +
  ggtitle("Chemistry (whole-transcriptome vs probe-based)")
print(p_condition | p_chem)

p_sample <- DimPlot(object, reduction = "umap.harmony", group.by = "sample_id",
                     shuffle = TRUE, seed = 1) +
  ggtitle("Sample identity (batch-correction check)") +
  theme(legend.text = element_text(size = 6)) +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))
print(p_sample)

message("==> Per-sample spatial cluster maps (10 samples) ...")
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
