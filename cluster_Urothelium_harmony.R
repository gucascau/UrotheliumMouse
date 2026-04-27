################################################################################
# Urothelium Harmony Reclustering
#
# Input:  integration_output/Urothelium_cells.rds  (265,354 cells, 3000 genes)
#         Seurat object produced by convert_Urothelium_to_rds.R;
#         the counts layer holds log-normalised X from the scVI h5ad.
#
# Steps:
#   1. Load Urothelium object
#   2. Set data layer = counts (already log-norm from scVI pipeline)
#   3. FindVariableFeatures (3000 HVGs)
#   4. ScaleData (regress pct_mt if present)
#   5. RunPCA (50 PCs)
#   6. RunHarmony (batch = sample_id + technology)
#   7. FindNeighbors + FindClusters
#   8. RunUMAP
#   9. Visualise → PDFs
#  10. Save Urothelium_harmony_integrated.rds
################################################################################

library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)
library(dplyr)

options(future.globals.maxSize = 4 * 1024^3)

log_mem <- function(label) {
  g  <- gc(verbose = FALSE)
  mb <- sum(g[, 2])
  message(sprintf("  [mem] %s: %.1f GB in use", label, mb / 1024))
}

# ── Paths ─────────────────────────────────────────────────────────────────────
DATA_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/"
OUT_DIR  <- file.path(DATA_DIR, "integration_output")
dir.create(OUT_DIR, showWarnings = FALSE)

setwd(paste0(OUT_DIR, "/UrotheliumAll/"))

UroDir <- paste0(OUT_DIR, "/UrotheliumAll/")
IN_PATH  <- file.path(OUT_DIR, "Urothelium_cells.rds")
OUT_PATH <- file.path(OUT_DIR, "Urothelium_harmony_integrated.rds")

# ── Parameters ────────────────────────────────────────────────────────────────
N_HVG        <- 3000
N_PCS        <- 20
HARMONY_VARS <- c("sample_id", "technology")
HARMONY_DIMS <- 1:20
RESOLUTION   <- 0.5


################################################################################
# STEP 1: Load Urothelium object
################################################################################

if (!file.exists(IN_PATH))
  stop("Input not found: ", IN_PATH)

message("Loading Urothelium_cells.rds ...")
so <- readRDS(IN_PATH)
message(sprintf("  Loaded: %d cells x %d genes", ncol(so), nrow(so)))
ASSAY <- DefaultAssay(so)
message(sprintf("  Default assay: %s", ASSAY))
log_mem("after load")
# so@meta.data %>% head()

# table(so@meta.data$technology)
# table(so@meta.data$condition)
# table(so@meta.data$gsm_id)
# table(so@meta.data$tissue)
# table(so@meta.data$paper)

# so@meta.data %>% filter(condition == "UUO_day2") %>% group_by(sample_id) %>% summarise(count = n())

# so@meta.data %>% filter(condition == "UUO_day2") %>% head()
# so@meta.data %>% filter(condition == "Sham") %>% head()
# so@meta.data %>% filter(condition == "Normal") %>% head()
# so@meta.data %>% filter(condition == "Development") %>% head()
# we need to seperate the KUDO urothelium sample to different sample Id
so@meta.data<- so@meta.data %>% mutate(
  sample_id = case_when(
    sample_id == "KudoUUOUrothelium" ~ Sample, 
    TRUE ~ sample_id),
  technology = case_when(
    Sample %in% c("GOF","LOF","Vehicle","Yoda") ~ "PIPseq", 
    Sample %in% c("Health","UUO_Day2","UUO_Day4","UUO_Day6","UUO_Day10","UUO_Day14") ~ "sci-RNA-seq3",
    Sample %in% c("GSM3827175_NMU_O_P","GSM3827176_NMU_O_D") ~ "10X",
    TRUE ~ technology),
  condition = case_when(
    Sample %in% c("GOF","LOF","Vehicle","Yoda") ~ "Urothelium_Organoid",
    Sample %in% c("Health") ~ "Healthy",
    Sample %in% c("UUO_Day2") ~ "UUO_2days",
    Sample %in% c("UUO_Day4") ~ "UUO_4days",
    Sample %in% c("UUO_Day6") ~ "UUO_6days",
    Sample %in% c("UUO_Day10") ~ "UUO_10days",
    Sample %in% c("UUO_Day14") ~ "UUO_14days",
    Sample %in% c("GSM3827175_NMU_O_P") ~ "Prolif_BlaUroOrganoid",
    Sample %in% c("GSM3827176_NMU_O_D") ~ "Dif_BlaUroOrganoid",
    condition %in% c("UUO_day2") ~ "UUO_2days",
    condition %in% c("UUO_day7") ~ "UUO_7days",
    condition %in% c("Sham") ~ "Healthy",
    condition %in% c("Normal") ~ "E18_5_Kidney",
    condition %in% c("Development") ~ "E9To13.5Gestation",
    TRUE ~ condition),
  paper = case_when(
    Sample %in% c("GOF","LOF","Vehicle","Yoda") ~ "Kudo2026",
    Sample %in% c("Health","UUO_Day2","UUO_Day4","UUO_Day6","UUO_Day10","UUO_Day14") ~ "PMID_36265491",
    Sample %in% c("GSM3827175_NMU_O_P","GSM3827176_NMU_O_D") ~ "PMID_31562298",
    TRUE ~ paper),
  tissue = case_when(
    Sample %in% c("GOF","LOF","Vehicle","Yoda") ~ "kidney",
    Sample %in% c("Health","UUO_Day2","UUO_Day4","UUO_Day6","UUO_Day10","UUO_Day14") ~ "kidney",
    Sample %in% c("GSM3827175_NMU_O_P") ~ "bladder",
    Sample %in% c("GSM3827176_NMU_O_D") ~ "bladder",
    TRUE ~ tissue),
  gsm_id = case_when(
    Sample %in% c("GOF","LOF","Vehicle","Yoda") ~ "CustomedKudoUrothelium",
    Sample %in% c("Health","UUO_Day2","UUO_Day4","UUO_Day6","UUO_Day10","UUO_Day14") ~ "GSE190887",
    Sample %in% c("GSM3827175_NMU_O_P") ~ "GSM3827175",
    Sample %in% c("GSM3827176_NMU_O_D") ~ "GSM3827176",
    TRUE ~ gsm_id)
)

table(so@meta.data$technology)
table(so@meta.data$condition)
table(so@meta.data$gsm_id)
table(so@meta.data$tissue)
table(so@meta.data$paper)


so@meta.data %>% filter(sample_id == "KudoUUOUrothelium") %>% group_by(species) %>% summarise(count = n())


table(so@meta.data$orig.ident)
DimPlot(so, group.by = "orig.ident", label = TRUE, repel = TRUE) + ggtitle("Urothelium coloured by orig.ident (pre-Harmony)")

################################################################################
# STEP 2: Set data layer = counts (X from h5ad is already log-normalised)
################################################################################

message("Setting data layer from counts (already log-norm) ...")
so[[ASSAY]]$data <- so[[ASSAY]]$counts
log_mem("after setting data layer")


################################################################################
# STEP 3: FindVariableFeatures
################################################################################

message(sprintf("FindVariableFeatures (nfeatures = %d) ...", N_HVG))
so <- FindVariableFeatures(so, selection.method = "vst",
                           nfeatures = N_HVG, verbose = FALSE)
message(sprintf("  HVGs selected: %d", length(VariableFeatures(so))))


################################################################################
# STEP 4: ScaleData
################################################################################

so[["pct_mt"]] <- PercentageFeatureSet(so, pattern = "^mt-")

so <- ScaleData(so,
                features        = VariableFeatures(so),
                vars.to.regress = "pct_mt",
                verbose         = FALSE)
message("  ScaleData done")
log_mem("after ScaleData")


################################################################################
# STEP 5: RunPCA
################################################################################

message(sprintf("RunPCA (%d PCs) ...", N_PCS))
so <- RunPCA(so, npcs = N_PCS, verbose = FALSE)
message("  PCA done")
log_mem("after PCA")

# # Free expression matrices before Harmony
# so[["RNA"]]$counts     <- NULL
# so[["RNA"]]$data       <- NULL
# so[["RNA"]]$scale.data <- NULL
# gc()
log_mem("after freeing expression data")


################################################################################
# STEP 6: RunHarmony
################################################################################

# Only use batch variables that actually exist in metadata
available_vars <- intersect(HARMONY_VARS, colnames(so@meta.data))
if (length(available_vars) == 0)
  stop("None of the requested Harmony variables found in metadata: ",
       paste(HARMONY_VARS, collapse = ", "))

message(sprintf("RunHarmony (batch = %s) ...",
                paste(available_vars, collapse = " + ")))

so <- RunHarmony(so, group.by.vars = HARMONY_VARS , reduction.save = "harmony",   
                 plot_convergence = TRUE, verbose = FALSE, dim.use = HARMONY_DIMS)
message("  Harmony done")
log_mem("after Harmony")


################################################################################
# STEP 7: FindNeighbors + FindClusters
################################################################################

message("FindNeighbors (annoy, k=20) ...")
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
# STEP 8: RunUMAP
################################################################################

message("RunUMAP on Harmony embedding ...")
so <- RunUMAP(so, reduction = "harmony", dims = HARMONY_DIMS,
              reduction.name = "umap_harmony", verbose = FALSE)
message("  UMAP done")


################################################################################
# STEP 9: Visualise
################################################################################

message("Generating UMAP plots ...")

make_plot <- function(grp, title, label = FALSE) {
  if (!grp %in% colnames(so@meta.data)) return(NULL)
  DimPlot(so, group.by = grp, reduction = "umap_harmony",
          label = label, repel = label, raster = TRUE) + ggtitle(title)
}

plots <- Filter(Negate(is.null), list(
  make_plot("seurat_clusters", "Clusters",   label = TRUE),
  make_plot("condition",       "By Condition"),
  make_plot("sample_id",       "By Sample"),
  make_plot("technology",      "By Technology"),
  make_plot("DataSet",         "By DataSet"),
  make_plot("paper",           "By Paper")
))

pdf(file.path(OUT_DIR, "Urothelium_UMAP_overview.pdf"), width = 18, height = 12)
print(wrap_plots(plots, ncol = 2))
dev.off()

# Marker gene feature plots
markers <- c("Krt5", "Krt14", "Krt20", "Krt8", "Krt18",
             "Upk2", "Upk3a", "Trp63", "Upk1a", "Upk1b")
present <- intersect(markers, rownames(so))
if (length(present) > 0) {
  # Restore data layer for FeaturePlot
  message("Restoring data layer for marker plots ...")
  so_tmp <- readRDS(IN_PATH)
  so[[ASSAY]]$data <- so_tmp[[ASSAY]]$counts
  rm(so_tmp); gc()

  pdf(file.path(OUT_DIR, "Urothelium_markers.pdf"), width = 18, height = 12)
  fp <- FeaturePlot(so, features = present, reduction = "umap_harmony",
                    ncol = 5, raster = TRUE)
  print(fp)
  dev.off()

  so[[ASSAY]]$data <- NULL
  gc()
}


################################################################################
# STEP 10: Save
################################################################################

message(sprintf("Saving to %s ...", OUT_PATH))
saveRDS(so, OUT_PATH)

message("\n===== Urothelium Harmony reclustering complete =====")
message(sprintf("Total cells: %s", format(ncol(so), big.mark = ",")))
message(sprintf("Clusters:    %d", length(unique(so$seurat_clusters))))
message(sprintf("Output:      %s", OUT_PATH))


################################################################################
# Briefly check the results
################################################################################

# read the saved object to confirm it loads correctly
so <- readRDS(OUT_PATH)

so@meta.data %>% head()


# we need to seperate the KUDO urothelium sample to different sample Id
so@meta.data$sample_id %>% table()
# filter the KudoUUOUrothelium
so@meta.data %>% filter(sample_id == "KudoUUOUrothelium") %>% group_by(Sample) %>% summarise(count = n())

# add the new sample id for so, when the sample is KudoUUOUrothelium, use the Sample column as the sample_id, otherwise keep the original sample_id



UrotheliumDimplotGroup<- DimPlot(so, group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
  ggtitle("Urothelium Harmony Reclustering: Clusters")

# save the plot as a PDF
ggsave(filename = file.path(UroDir, "Urothelium_UMAP_clusters_ByGroup.pdf"), plot = UrotheliumDimplotGroup, width = 8, height = 6)

UrotheliumDimplotGroupSampleID<- DimPlot(so, group.by = "sample_id", label = TRUE, repel = TRUE) +
  ggtitle("Urothelium Harmony Reclustering: Clusters")

# save the plot as a PDF
ggsave(filename = file.path(UroDir, "Urothelium_UMAP_GroupbySampleID.pdf"), plot = UrotheliumDimplotGroupSampleID, width = 8, height = 6)


UrotheliumDimplotGrouptissueSplitSample<-DimPlot(so, group.by = "tissue", split.by = "sample_id",label = TRUE, repel = TRUE, ncol = 4) +
  ggtitle("Urothelium Harmony Reclustering: tissue - split by sample_id")

# save the plot as a PDF
ggsave(filename = file.path(UroDir, "UrotheliumDimplotGrouptissueSplitSample.pdf"), plot = UrotheliumDimplotGrouptissueSplitSample, width = 10, height = 10)

UrotheliumDimplotGroupSplitSample<-DimPlot(so, group.by = "seurat_clusters", split.by = "sample_id",label = TRUE, repel = TRUE, ncol = 4) +
  ggtitle("Urothelium Harmony Reclustering: Clusters- split by sample_id")

# save the plot as a PDF
ggsave(filename = file.path(UroDir, "Urothelium_UMAP_clusters_ByGroup_Splitbysample.pdf"), plot = UrotheliumDimplotGroupSplitSample, width = 20, height = 20)


# The samples from the development dataset are quite complex, we will ignore these developmental datasets
so1 <- subset(so, subset = condition != "E9To13.5Gestation" & condition != "E18_5_Kidney")

# Generate the UMAP plot coloured by conditions, 
message("Generating UMAP plots for dataset without developmental...")

make_plot_dimplot <- function(obj, grp, title, label = FALSE) {
  if (!grp %in% colnames(obj@meta.data)) return(NULL)
  DimPlot(obj, group.by = grp, reduction = "umap_harmony",
          label = label, repel = label, raster = TRUE) + ggtitle(title)
}

plots <- Filter(Negate(is.null), list(
  make_plot_dimplot(so1, "seurat_clusters", "Clusters",   label = TRUE),
  make_plot_dimplot(so1, "condition",       "By Condition"),
  make_plot_dimplot(so1, "sample_id",       "By Sample"),
  make_plot_dimplot(so1, "technology",      "By Technology"),
  make_plot_dimplot(so1, "DataSet",         "By DataSet"),
  make_plot_dimplot(so1, "paper",           "By Paper"),
  make_plot_dimplot(so1, "tissue",           "By tissue")
))

pdf(file.path(OUT_DIR, "Urothelium_nodevelopment_UMAP_overview.pdf"), width = 18, height = 12)
print(wrap_plots(plots, ncol = 2))
dev.off()


# Generate the UMAP plot coloured by conditions, split by condition, without developmental datasets
message("Generating UMAP plots for dataset without developmental, split by the conditions...")

make_plot_dimplotsplit <- function(obj, grp, title, label = FALSE) {
  if (!grp %in% colnames(obj@meta.data)) return(NULL)
  DimPlot(obj, group.by = grp, split.by = grp, reduction = "umap_harmony",
          label = label, repel = label, raster = TRUE) + ggtitle(title)
}

plots <- Filter(Negate(is.null), list(
  make_plot_dimplotsplit(so1, "seurat_clusters", "Clusters",   label = TRUE),
  make_plot_dimplotsplit(so1, "condition",       "Split By Condition"),
  make_plot_dimplotsplit(so1, "sample_id",       "Split By Sample"),
  make_plot_dimplotsplit(so1, "technology",      "Split By Technology"),
  make_plot_dimplotsplit(so1, "paper",           "Split By Paper"),
  make_plot_dimplotsplit(so1, "tissue",           "Split By tissue")
))

# generate a loop file to save the plot
# save the file based on the number of plots. 
# # calculate the number of plots based on the factors, including condition, sample_id, technology, paper, tissue.
# so1@meta.data$condition %>% unique() %>% length()
# so1@meta.data$Sample %>% unique() %>% length()
# so1@meta.data$technology %>% unique() %>% length()
# so1@meta.data$paper %>% unique() %>% length()


# for (i in 1:length(plots)) {
#   plot <- plots[[i]]
#   title <- plot$labels$title
#   filename <- paste0("Urothelium_nodevelopment_UMAP", gsub(" ", "_", title), ".pdf")
#   ggsave(filename = file.path(OUT_DIR, filename), plot = plot, width = 20, height = 5)
# }


# Map each plot title to its metadata column for counting levels
group_cols <- c(
  "By Condition"  = "condition",
  "By Sample"     = "Sample",
  "By Technology" = "technology",
  "By Paper"      = "paper",
  "By Tissue"     = "tissue"
)

for (i in seq_along(plots)) {
  plot  <- plots[[i]]
  title <- plot$labels$title
  col   <- group_cols[title]

  n_groups <- if (!is.na(col) && col %in% colnames(so1@meta.data)) {
    length(unique(so1@meta.data[[col]]))
  } else 10L

  # UMAP body ~7in; legend ~0.25in per item, capped at 2 legend columns
  width  <- 10 + min(n_groups, 30) * 0.25 + if (n_groups > 30) 2 else 0
  height <- max(7, ceiling(n_groups / 2) * 0.2 + 6)

  filename <- paste0("Urothelium_nodevelopment_UMAP_", gsub(" ", "_", title), ".pdf")
  ggsave(file.path(OUT_DIR, filename), plot = plot,
         width = width, height = height)
}

# plot the biomarkers that without development
markers <- c("Krt5", "Krt14", "Krt20", "Krt8", "Krt18",
             "Upk2", "Upk3a",  "Upk1a", "Upk1b","Trp63")
if (length(markers) > 0) {
  # Restore data layer for FeaturePlot

  pdf(file.path(OUT_DIR, "Urothelium_markers_nodevelopment.pdf"), width = 18, height = 12)
  fp <- FeaturePlot(so1, features = markers, reduction = "umap_harmony",
                    ncol = 5, raster = TRUE)
  print(fp)
  dev.off()
}


# generate the dotplot across various metadata groups for the marker genes, without developmental datasets
metadata_groups <- c("condition", "Sample", "technology", "paper", "tissue")
DotPlot(so1, features = markers, group.by = "condition") + RotatedAxis() + ggtitle("DotPlot of markers by condition")
