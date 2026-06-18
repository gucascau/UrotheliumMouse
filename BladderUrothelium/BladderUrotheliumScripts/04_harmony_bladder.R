################################################################################
# Script 04: Seurat + Harmony Integration — Bladder Urothelium
#
# Steps:
#   1. Load all *_qc.rds from seurat_objects/
#   2. Merge into one Seurat object
#   3. NormalizeData + FindVariableFeatures (3000 HVGs, batch-aware)
#   4. JoinLayers + ScaleData (regress pct_mt) + RunPCA (50 PCs)
#   5. RunHarmony (batch = sample_id + technology)
#   6. FindNeighbors + FindClusters
#   7. RunUMAP
#   8. Visualize → PDFs
#   9. Save bladder_harmony_integrated.rds
#
# Requires: *_qc.rds produced by 02_qc_bladder.R
################################################################################

library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)

options(future.globals.maxSize = 8 * 1024^3)

log_mem <- function(label) {
  g  <- gc(verbose = FALSE)
  mb <- sum(g[, 2])
  message(sprintf("  [mem] %s: %.1f GB in use", label, mb / 1024))
}

DATA_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/BladderUrothelium"
OBJ_DIR  <- file.path(DATA_DIR, "seurat_objects")
OUT_DIR  <- file.path(DATA_DIR, "integration_output")
dir.create(OUT_DIR, showWarnings = FALSE)

OUT_PATH <- file.path(OUT_DIR, "bladder_harmony_integrated.rds")

N_HVG        <- 3000
N_PCS        <- 50
HARMONY_VARS <- c("sample_id", "technology")
HARMONY_DIMS <- 1:30
RESOLUTION   <- 0.5


################################################################################
# STEP 1: Load all QC'd objects
################################################################################

qc_files   <- sort(list.files(OBJ_DIR, pattern = "_qc\\.rds$",
                               full.names = TRUE))
sample_ids <- sub("_qc\\.rds$", "", basename(qc_files))
n          <- length(qc_files)
message(sprintf("Loading %d QC'd Seurat objects...", n))

seurat_list <- lapply(seq_along(qc_files), function(i) {
  message(sprintf("  [%d/%d] %s", i, n, sample_ids[i]))
  so <- readRDS(qc_files[i])
  # Ensure RNA is active assay; drop SCT from DoubletFinder preprocessing
  if (!"RNA" %in% Assays(so)) {
    so <- RenameAssays(so, assay.name = DefaultAssay(so),
                       new.assay.name = "RNA")
  }
  DefaultAssay(so) <- "RNA"
  for (a in setdiff(Assays(so), "RNA")) so[[a]] <- NULL
  so <- JoinLayers(so)
  so
})
names(seurat_list) <- sample_ids
log_mem("after loading all objects")


################################################################################
# STEP 2: Merge
################################################################################

message("Merging objects...")
merged <- merge(seurat_list[[1]],
                y          = seurat_list[-1],
                add.cell.ids = sample_ids,
                merge.data   = FALSE)
rm(seurat_list); gc()
message(sprintf("  Merged: %d cells × %d genes", ncol(merged), nrow(merged)))
log_mem("after merge")


################################################################################
# STEP 3: Normalize + batch-aware HVG selection
################################################################################

message("NormalizeData...")
merged <- NormalizeData(merged, normalization.method = "LogNormalize",
                        scale.factor = 1e4, verbose = FALSE)

message(sprintf("FindVariableFeatures (nfeatures = %d, batch-aware)...", N_HVG))
merged <- FindVariableFeatures(merged,
                               selection.method = "vst",
                               nfeatures        = N_HVG,
                               verbose          = FALSE)
message(sprintf("  HVGs selected: %d", length(VariableFeatures(merged))))

# Compute mitochondrial percentage if not already present
if (!"pct_mt" %in% colnames(merged@meta.data)) {
  merged[["pct_mt"]] <- PercentageFeatureSet(merged, pattern = "^mt-|^MT-")
}


################################################################################
# STEP 4: JoinLayers + ScaleData + RunPCA
################################################################################

message("JoinLayers...")
merged <- JoinLayers(merged)
log_mem("after JoinLayers")

message("ScaleData (regress pct_mt, HVGs only)...")
merged <- ScaleData(merged,
                    features        = VariableFeatures(merged),
                    vars.to.regress = "pct_mt",
                    verbose         = FALSE)
log_mem("after ScaleData")

message(sprintf("RunPCA (%d PCs)...", N_PCS))
merged <- RunPCA(merged, npcs = N_PCS, verbose = FALSE)
log_mem("after PCA")

################################################################################
# STEP 5: RunHarmony
################################################################################

available_vars <- intersect(HARMONY_VARS, colnames(merged@meta.data))
if (length(available_vars) == 0)
  stop("No Harmony batch variables found in metadata: ",
       paste(HARMONY_VARS, collapse = ", "))

message(sprintf("RunHarmony (batch = %s)...",
                paste(available_vars, collapse = " + ")))
merged <- RunHarmony(merged,
                     group.by.vars  = available_vars,
                     reduction      = "pca",
                     dims.use       = HARMONY_DIMS,
                     reduction.save = "harmony",
                     plot_convergence = TRUE,
                     verbose        = FALSE)
message("  Harmony done")
log_mem("after Harmony")

# Free scale.data after Harmony (harmony needs it internally with Seurat v5)
merged[["RNA"]]$scale.data <- NULL
gc()
log_mem("after freeing scale.data")


################################################################################
# STEP 6: FindNeighbors + FindClusters
################################################################################

message("FindNeighbors (annoy, k = 20)...")
merged <- FindNeighbors(merged,
                        reduction    = "harmony",
                        dims         = HARMONY_DIMS,
                        nn.method    = "annoy",
                        k.param      = 20,
                        annoy.metric = "euclidean",
                        n.trees      = 50,
                        verbose      = FALSE)
log_mem("after FindNeighbors")

message(sprintf("FindClusters (resolution = %.2f)...", RESOLUTION))
merged <- FindClusters(merged, resolution = RESOLUTION, verbose = FALSE)
message(sprintf("  Clusters: %d", length(unique(merged$seurat_clusters))))
gc()


################################################################################
# STEP 7: RunUMAP
################################################################################

message("RunUMAP on Harmony embedding...")
merged <- RunUMAP(merged,
                  reduction      = "harmony",
                  dims           = HARMONY_DIMS,
                  reduction.name = "umap_harmony",
                  verbose        = FALSE)
message("  UMAP done")


################################################################################
# STEP 8: Visualize
################################################################################

message("Generating UMAP plots...")

make_plot <- function(grp, title, label = FALSE) {
  if (!grp %in% colnames(merged@meta.data)) return(NULL)
  DimPlot(merged, group.by = grp, reduction = "umap_harmony",
          label = label, repel = label, raster = TRUE) + ggtitle(title)
}

plots <- Filter(Negate(is.null), list(
  make_plot("seurat_clusters", "Clusters",        label = TRUE),
  make_plot("condition",       "By Condition"),
  make_plot("sample_id",       "By Sample"),
  make_plot("technology",      "By Technology"),
  make_plot("gsm_id",          "By GSM ID"),
  make_plot("paper",           "By Paper")
))

pdf(file.path(OUT_DIR, "bladder_harmony_UMAP_overview.pdf"),
    width = 18, height = 12)
print(wrap_plots(plots, ncol = 2))
dev.off()

# Urothelial marker feature plots
markers <- c("Krt5", "Krt14", "Krt20", "Krt8", "Krt18",
             "Upk2", "Upk3a", "Trp63", "Upk1a", "Upk1b",
             "Cd44", "Vim", "Acta2", "Pecam1", "Ptprc")
present <- intersect(markers, rownames(merged))
if (length(present) > 0) {
  pdf(file.path(OUT_DIR, "bladder_harmony_markers.pdf"),
      width = 18, height = ceiling(length(present) / 5) * 4)
  print(FeaturePlot(merged, features = present, reduction = "umap_harmony",
                    ncol = 5, raster = TRUE))
  dev.off()
}


################################################################################
# STEP 9: Save
################################################################################

message(sprintf("Saving to %s ...", OUT_PATH))
saveRDS(merged, OUT_PATH)

message("\n===== Bladder Harmony integration complete =====")
message(sprintf("Total cells: %s", format(ncol(merged), big.mark = ",")))
message(sprintf("Clusters:    %d", length(unique(merged$seurat_clusters))))
message(sprintf("Samples:     %d", length(unique(merged$sample_id))))
message(sprintf("Output:      %s", OUT_PATH))
