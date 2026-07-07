library(Seurat)
library(ggplot2)

OUT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/SpatialScripts/output"
region_obj <- readRDS(file.path(OUT_DIR, "VisiumHD_region_StromalUrothelialMyeloid.rds"))
message("Assays: ", paste(Assays(region_obj), collapse=", "))
message("Images: ", paste(Images(region_obj), collapse=", "))
message("dims: ", paste(dim(region_obj), collapse=" x "))

gene <- "Krt5"
image_name <- "slice1.008um.2"

# Fix A: SpatialFeaturePlot with explicit alpha range so zero-expression bins fade out
pA <- SpatialFeaturePlot(region_obj,
  features = gene, images = image_name, crop = TRUE,
  pt.size.factor = 3, image.alpha = 1,
  alpha = c(0, 1)
) + ggtitle(paste("Fix A: alpha=c(0,1) -", gene)) +
  scale_fill_gradientn(colors = c("transparent","#FFFFD4","#FED98E","#FE8929","#CC4C02"), na.value="transparent")

ggsave(file.path(OUT_DIR, "test_fixA_alpha.png"), pA, width=6, height=5, dpi=150)

# Fix B: ImageFeaturePlot as suggested by user
pB <- tryCatch({
  ImageFeaturePlot(region_obj,
    features = gene, fov = image_name,
    dark.background = FALSE,
    cols = c("lightgrey", "red")
  ) + ggtitle(paste("Fix B: ImageFeaturePlot -", gene))
}, error = function(e) { message("ImageFeaturePlot error: ", conditionMessage(e)); NULL })

if (!is.null(pB)) {
  ggsave(file.path(OUT_DIR, "test_fixB_imagefeatureplot.png"), pB, width=6, height=5, dpi=150)
} else {
  message("Fix B failed, skipping save")
}

message("Done.")
