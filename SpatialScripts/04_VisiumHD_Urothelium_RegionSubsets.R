################################################################################
# 04_VisiumHD_Urothelium_RegionSubsets.R
#
# Subset the RCTD-deconvolved VisiumHD object (both kidneys) into three
# manually-selected regions, save each as its own RDS, and plot urothelium
# biomarker spatial expression for each region.
#
# Region CSVs and their verified sample of origin (see
# resolve_3side_urothelium_region.py and inline notes below for how this was
# confirmed empirically — the barcode formats/sample assignment are NOT what
# an earlier half-finished exploration block in
# 01_VisiumHD_integrate_harmony_sketch.R assumed):
#
#   stromal–urothelial–myeloid .csv          -> kidney5p (slice1.008um.2)
#   Selected_UrotheliumRegionV2.csv           -> kidney5p (slice1.008um.2)
#   Selected_3SideUrotheliumRegion_V2.csv     -> kidney3p (slice1.008um)
#     (segmented-cell IDs, pre-resolved to square_008um bins via
#      Selected_3SideUrotheliumRegion_V2_resolved_008um_bins.txt)
#
# Cell name prefixes in the integrated object (from the merge() call in
# 01_VisiumHD_integrate_harmony_sketch.R, add.cell.ids = c("kidney3p","kidney")):
#   kidney3p sample -> "kidney3p_" prefix, image "slice1.008um"
#   kidney5p sample -> "kidney_"    prefix, image "slice1.008um.2"
#
# Input:  VisiumHD_kidney_both_deconvolved.rds
#         stromal–urothelial–myeloid .csv
#         Selected_UrotheliumRegionV2.csv
#         Selected_3SideUrotheliumRegion_V2_resolved_008um_bins.txt
# Output: VisiumHD_region_StromalUrothelialMyeloid.rds (+ PDF)
#         VisiumHD_region_UrotheliumRegionV2.rds (+ PDF)
#         VisiumHD_region_3SideUrothelium.rds (+ PDF)
################################################################################

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
# ── Paths ─────────────────────────────────────────────────────────────────────
HD_BASE <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSpatialData/VisiumHD"
OUT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/SpatialScripts/output"

OBJECT_RDS <- file.path(OUT_DIR, "VisiumHD_kidney_both_deconvolved.rds")

# Urothelial marker panel used throughout this project (see UrothelialMarkers
# in 01_VisiumHD_integrate_harmony_sketch.R line ~669)


RequestedCategoriesName <- "UrotheliumMarkers"
Requested_MARKERS <- c(
  "Krt5", "Krt14", "Krt20", "Trp63",
  "Upk1a", "Upk1b", "Upk2", "Upk3a", "Upk3b",
  "Foxa1", "Gata3", "Pparg",
  "Krt8", "Krt18", "Krt19"
)

# RequestedCategoriesName <- "Spp1Genes"
# Requested_MARKERS<- c("Spp1","Cd44","Itgav","Itgb1","Itgb3","Itgb5","Itgb6","Itgb8")

RequestedCategoriesName <- "RenalInjuryMarkers"
Requested_MARKERS <- c(
  "Spp1","Havcr1", "Lcn2", "Clu", "Cd44","Fn1", "Vcam1", "Icam1"
)


# change the legend title to "PerExp" and "AveExp" for the dotplot
# NeighborCells markers: 
RequestedCategoriesName <- "NeighborCellsMarkers"

Requested_MARKERS <- c("Aqp2","Aqp3","Col1a2","Col3a1","Acta2","Tagln","Cdh5","Pecam1","Cd3e", "Cd4", "Cd8a", "Cd19", "Cd79a", "C1qa","Itgam", "Adgre1", "Ly6g")


RequestedCategoriesName <- "EnrichedMarkers"

Requested_MARKERS <- c("Gsdmc2","Fxyd3","Snx31","Sprr1a","Foxq1","Psca","S100a6 
","Lgals3")

RequestedCategoriesName <- "UpperRenalUroGenes"

Requested_MARKERS <- c("Pax8","Pax2","Glis3","Fgfr1","Pkhd1","Bicc1","Magi1","Cgnl1","Ptpn14", "Col4a3", "Col4a4", "Col4a5")


# ── Load deconvolved object ────────────────────────────────────────────────────
message("==> Loading VisiumHD deconvolved object ...")
object <- readRDS(OBJECT_RDS)
# Object also carries a "sketch" assay (5000-bin/sample LeverageScore subsample
# used for fast Harmony clustering in 01_VisiumHD_integrate_harmony_sketch.R).
# DefaultAssay affects Cells(object)/subset() below — must be the full-resolution
# assay, or barcode matching silently checks against only the 10,000-bin sketch
# instead of the ~925,000 full-resolution bins (this was the actual cause of the
# "barcodes requested: 9957 | matched: 179"-type mismatches).
DefaultAssay(object) <- "Spatial.008um"
all_images <- Images(object)
message("Image slots found: ", paste(all_images, collapse = ", "))
# before doing anything I want to change the


# Generate the cell type within nich 9
object @meta.data %>% filter (seurat_cluster.harmony.projected == "9") %>% head()
# we extract the RCTD columns
rctd_cols<- grep("^rctd_", colnames(object@meta.data), value = T)[1:36]
object @meta.data %>% select(all_of(rctd_cols)) 

# I want to combine the rctd_Asc-Vasa-Recta, rctd_Desc-Vasa-Recta, rctd_Vas-Efferens, rctd_Vas-Afferens to rctd_endo, change the name of rctd_Per -> rctd_Peri/VSMC
# ── Collapse related RCTD cell types for interpretability ─────────────────────
VASCULAR_SUBTYPES <- c("rctd_Asc-Vasa-Recta", "rctd_Desc-Vasa-Recta",
                        "rctd_Vas-Efferens", "rctd_Vas-Afferens")

collapse_rctd_celltypes <- function(df) {
  vascular_present <- intersect(VASCULAR_SUBTYPES, colnames(df))
  if (length(vascular_present) > 0) {
    # No na.rm: a bin's rctd_* weights are either all NA (RCTD-excluded bin)
    # or all numeric together, so plain rowSums preserves that NA-ness
    # instead of silently turning excluded bins into a fabricated 0.
    df$rctd_Endo <- rowSums(df[, c("rctd_Endo", vascular_present), drop = FALSE])
    df <- df[, !(colnames(df) %in% vascular_present)]
  }
  if ("rctd_Per" %in% colnames(df)) {
    colnames(df)[colnames(df) == "rctd_Per"] <- "rctd_Peri/VSMC"
  }
  df
}
meta <- object@meta.data
meta <- collapse_rctd_celltypes(meta)

rctd_cols <- grep("^rctd_", colnames(meta), value = TRUE)
rctd_cols <- setdiff(rctd_cols, "rctd_dominant_celltype")

# Then I added this into the object
object@meta.data %>% head()
meta %>% head()
# Normalize weights so each spot's cell-type proportions sum to 1

# Subset to the cluster of interest
cluster_id <- "9"  # change to your cluster of interest
 # cluster_weights <- # Subset metadata to cells in the cluster of interest
composition_df <- meta %>%
  filter(seurat_cluster.harmony.projected == cluster_id) %>%
  select(sample_id, all_of(rctd_cols)) %>%
  group_by(sample_id) %>%
  summarise(across(all_of(rctd_cols), \(x) mean(x, na.rm = TRUE))) %>%
  pivot_longer(cols = all_of(rctd_cols), names_to = "cell_type", values_to = "proportion") %>%
  mutate(cell_type = sub("^rctd_", "", cell_type))
## Optional: collapse rare types into "Other" per sample
composition_df <- composition_df %>%
  mutate(cell_type = ifelse(proportion < 0.01, "Other", cell_type)) %>%
  group_by(sample_id, cell_type) %>%
  summarise(proportion = sum(proportion), .groups = "drop")

# 
ggplot(composition_df, aes(x = sample_id, y = proportion, fill = cell_type)) +
  geom_col(position = "stack") +
  labs(title = paste("Cell type composition by sample — Cluster", cluster_id),
       x = NULL, y = "Mean proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Or as a pie chart
Nich9CelltypePieCharDonot<- ggplot(composition_df, aes(x = 2, y = proportion, fill = cell_type)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +  # creates the donut hole
  facet_wrap(~ sample_id) +
  theme_void() 
  #+
  # labs(title = paste("Cluster", cluster_id, "— Cell type composition by sample"))

# save this pie chart for the nich9
ggsave(file.path(OUT_DIR, "XeniumHD_Nich9CelltypePieCharDonot.pdf"), plot= Nich9CelltypePieCharDonot, height =4, width =7)

ggplot(composition_df, aes(x = "", y = proportion, fill = cell_type)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  facet_wrap(~ sample_id) +
  theme_void() +
  labs(title = paste("Cluster", cluster_id, "— Cell type composition by sample")) +
  theme(strip.text = element_text(size = 11, face = "bold"))

markers_present <- intersect(Requested_MARKERS, rownames(object))

# RequestedCategoriesName <- "Spp1Genes"
# Spp1Genes<- c("Spp1","Cd44","Itgav","Itgb1","Itgb3","Itgb5","Itgb6","Itgb8")
# markers_present <- intersect(Spp1Genes, rownames(object))

message(sprintf("Urothelial markers present in object: %d / %d",
                length(markers_present), length(Requested_MARKERS)))
message("  ", paste(markers_present, collapse = ", "))

# ── Helper: subset by barcode list, save RDS, plot markers ────────────────────
subset_region <- function(barcodes, prefix, image_name, region_name) {
  message(sprintf("\n==> Region: %s", region_name))

  cell_ids <- paste0(prefix, barcodes)
  matched  <- intersect(cell_ids, Cells(object))
  message(sprintf("  Barcodes requested: %d | matched in object: %d",
                  length(cell_ids), length(matched)))

  if (length(matched) == 0) {
    stop(sprintf(
      "No cells matched for region '%s' (prefix '%s', image '%s'). ",
      region_name, prefix, image_name),
      "Barcode prefix or sample assignment is likely wrong — do not proceed silently."
    )
  }

  region_obj <- subset(object, cells = matched)
  region_obj <- JoinLayers(region_obj)

  out_rds <- file.path(OUT_DIR, sprintf("VisiumHD_region_%s_%s.rds", region_name, RequestedCategoriesName))
  saveRDS(region_obj, out_rds)
  message("  Saved: ", out_rds)

  out_pdf <- file.path(OUT_DIR, sprintf("VisiumHD_region_%s_%s.pdf", region_name, RequestedCategoriesName))
  pdf(out_pdf, width = 3.5, height = 4)
  for (gene in markers_present) {
    p <- SpatialFeaturePlot(region_obj,
      features       = gene,
      images         = image_name,
      crop           = TRUE,
      pt.size.factor = 12,
      image.alpha    = 1,
      alpha = c(0.1, 1)
    ) + ggtitle(paste(region_name, "-", gene)) +
      scale_fill_gradientn(
        colors   = c("transparent", "#FFFFD4", "#FED98E", "#FE8929", "#CC4C02"),
        na.value = "transparent"
      )
    print(p)
  }
  dev.off()

  message("  Saved: ", out_pdf)

  region_obj
}

# ── Region 1: the top urothelium region, that renal pelviscalyceal urotheliual line (kidney5p) ───────────────────────────
csv1 <- read.csv(file.path(HD_BASE, "Kidney1TopUroSection.csv"))
region1 <- subset_region(
  barcodes    = csv1$Barcode,
  prefix      = "kidney_",
  image_name  = "slice1.008um.2",
  region_name = "Pelviscalyceal_urothelium"
)
rm(region1); gc()

# ── Region 2: the middle urothelium region, that rrenal pelvis calyx lumen line (kidney5p) ────────────────────────
csv2 <- read.csv(file.path(HD_BASE, "Kidney1MiddleUroSection.csv"))
region2 <- subset_region(
  barcodes    = csv2$Barcode,
  prefix      = "kidney_",
  image_name  = "slice1.008um.2",
  region_name = "PelvisCalyxLumen_urothelium"
)
rm(region2); gc()


# ── Region 3: the bottom urothelium region, Suburothelial peri-urothelial region (kidney5p) ────────────────────────
csv3 <- read.csv(file.path(HD_BASE, "Kidney1BottomUroSection.csv"))
region3 <- subset_region(
  barcodes    = csv3$Barcode,
  prefix      = "kidney_",
  image_name  = "slice1.008um.2",
  region_name = "Periurothelial_urothelium"
)
rm(region3); gc()


# ── Region 4: 3-side middleurothelium region (kidney3p, resolved from cell IDs) ─────
resolved_bins_file <- file.path(HD_BASE, "Kidney3P_LeftUrothelium.csv")

bins4 <- read.csv(resolved_bins_file)
region4 <- subset_region(
  barcodes    = bins4$Barcode,
  prefix      = "kidney3p_",
  image_name  = "slice1.008um",
  region_name = "3SideLeftUrothelium"
)
rm(region4); gc()

# ── Region 5: 3-side right urothelium region (kidney3p, resolved from cell IDs) ─────
resolved_bins_file <- file.path(HD_BASE, "Kidney3P_RightUrothelium.csv")
bins5 <- read.csv(resolved_bins_file)
region5 <- subset_region(
  barcodes    = bins5$Barcode,
  prefix      = "kidney3p_",
  image_name  = "slice1.008um",
  region_name = "3SideRightUrothelium"
)
rm(region5); gc()

message("\n==> Done.")


# We incoporated all these five reigons into one Seurat object for downstream analysis

Pelviscalyceal_urothelium <- readRDS(file.path(OUT_DIR, "VisiumHD_region_Pelviscalyceal_urothelium_UrotheliumMarkers.rds"))
PelvisCalyxLumen_urothelium <- readRDS(file.path(OUT_DIR, "VisiumHD_region_PelvisCalyxLumen_urothelium_UrotheliumMarkers.rds"))
Periurothelial_urothelium <- readRDS(file.path(OUT_DIR, "VisiumHD_region_Periurothelial_urothelium_UrotheliumMarkers.rds"))
ThreeSideLeftUrothelium <- readRDS(file.path(OUT_DIR, "VisiumHD_region_3SideLeftUrothelium_UrotheliumMarkers.rds"))
ThreeSideRightUrothelium <- readRDS(file.path(OUT_DIR, "VisiumHD_region_3SideRightUrothelium_UrotheliumMarkers.rds"))

# we can merge all these five regions into one Seurat object for downstream analysis
UrotheliumSelectedRegions <- merge(Pelviscalyceal_urothelium, y = c(PelvisCalyxLumen_urothelium, Periurothelial_urothelium, ThreeSideLeftUrothelium, ThreeSideRightUrothelium), add.cell.ids = c("Pelviscalyceal", "PelvisCalyxLumen", "Periurothelial", "ThreeSideLeft", "ThreeSideRight"), project = "UrotheliumSelectedRegions")

# add the metadata to the merged Seurat object
UrotheliumSelectedRegions$Region <- sapply(strsplit(colnames(UrotheliumSelectedRegions), "_"), `[`, 1)

# save the merged Seurat object
saveRDS(UrotheliumSelectedRegions, file.path(OUT_DIR, "VisiumHD_UrotheliumSelectedFiveRegions_UrotheliumMarkers.rds"))

UrotheliumSelectedRegions <-readRDS(file.path(OUT_DIR, "VisiumHD_UrotheliumSelectedFiveRegions_UrotheliumMarkers.rds"))
# We only concentrate on the cluster 9 for these integrated Seurat object, which is the urothelium cluster. We will subset the cluster 9 and save it as a new Seurat object for downstream analysis.
OnlyUrotheliumCluster <- subset(UrotheliumSelectedRegions, idents = "9")

# change the legend title to "PerExp" and "AveExp" for the dotplot
RequestedGeneDotPlotInUro <-DotPlot(OnlyUrotheliumCluster, features = Requested_MARKERS, group.by = "Region", cols = c("lightgrey", "red"), dot.scale = 8) +  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(
    color = guide_colorbar(title = "AveExp"),
    size  = guide_legend(title = "PerExp")
  )

# save the dotplot as a pdf
ggsave(file.path(OUT_DIR, sprintf("DotPlot_UrotheliumSelectedFiveRegions_%s.pdf", RequestedCategoriesName)), RequestedGeneDotPlotInUro, width = 6, height = 3.5)

OnlyUrotheliumCluster<- JoinLayers(OnlyUrotheliumCluster)
Idents(OnlyUrotheliumCluster) <- OnlyUrotheliumCluster$Region
# Identify the top 5 markers for each region in the merged Seurat object
RegionalDEGs <- FindAllMarkers(OnlyUrotheliumCluster, only.positive = TRUE, min.pct = 0.2, logfc.threshold = 0.2, test.use = "MAST") 

# Draw the heatmap for the top 5 markers for each region in the merged Seurat object
Top5MarkersPerRegion <- RegionalDEGs %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
Top5MarkersPerRegionHeatmap <- DoHeatmap(OnlyUrotheliumCluster, features = Top5MarkersPerRegion$gene, group.by = "Region", size = 3) + scale_fill_gradientn(colors = c("blue", "white", "red"), na.value = "white") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# save the heatmap as a pdf
ggsave(file.path(OUT_DIR, sprintf("Heatmap_Top5MarkersPerRegion_UrotheliumSelectedFiveRegions_%s.pdf", RequestedCategoriesName)), Top5MarkersPerRegionHeatmap, width = 6, height = 4)

