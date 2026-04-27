################################################################################
# 05b_harmony_recluster_renal_epithelial.R
#
# Harmony reclustering of the full renal epithelial subset (tubule, urothelium,
# collecting duct, glomerular epithelial).
#
# Input : output/RenalUrothelium_renal_epithelial.rds
#         counts layer = log-normalised X (same convention as 05_harmony_recluster.R)
#
# Steps:
#   1.  Load Seurat object
#   2.  Set data layer = counts (already log-norm)
#   3.  Compute pct_mt + FindVariableFeatures (3 000 HVGs)
#   4.  ScaleData (regress pct_mt)
#   5.  RunPCA (30 PCs — larger than urothelium-only pipeline to capture
#               the extra biological diversity across all four compartments)
#   6.  RunHarmony (sample_id; +technology/source when >1 level)
#   7.  FindNeighbors + FindClusters (res = 0.5)
#   8.  RunUMAP
#   9.  UMAP + marker FeaturePlots → PDFs
#  10.  Save output/RenalUrothelium_renal_epithelial_harmony.rds
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

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR  <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
OUT_DIR   <- file.path(BASE_DIR, "RenalUrotheliumScripts", "output")
PLOT_DIR  <- file.path(BASE_DIR, "RenalUrotheliumScripts", "plots", "renal_epithelial")
dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

IN_PATH  <- file.path(OUT_DIR, "RenalUrothelium_renal_epithelial.rds")
OUT_PATH <- file.path(OUT_DIR, "RenalUrothelium_renal_epithelial_harmony.rds")

# ── Parameters ────────────────────────────────────────────────────────────────
N_HVG        <- 3000
N_PCS        <- 30       # more PCs than urothelium-only to capture 4-compartment diversity
HARMONY_DIMS <- 1:30
RESOLUTION   <- 0.5

# ── Per-compartment marker genes for feature plots ────────────────────────────
MARKERS <- list(
  Tubule = c(
    "Lrp2", "Slc34a1", "Cubn", "Hnf4a",          # proximal tubule
    "Umod", "Slc12a1", "Cldn16",                   # loop of Henle / TAL
    "Slc12a3", "Pvalb",                            # distal convoluted tubule
    "Calb1", "Six2"                                # CNT / nephron progenitor
  ),
  Urothelium = c(
    "Upk1a", "Upk1b", "Upk2", "Upk3a",
    "Krt5", "Krt14", "Krt20", "Krt8", "Krt18", "Trp63"
  ),
  Collecting_duct = c(
    "Aqp2", "Avpr2", "Hsd11b2",                   # principal cells
    "Atp6v1b1", "Foxi1", "Slc4a1", "Slc26a4"      # intercalated cells
  ),
  Glomerular_epithelial = c(
    "Nphs1", "Nphs2", "Podxl", "Wt1", "Synpo",    # podocytes
    "Pax8", "Cldn1", "Akr1b7"                      # parietal epithelial cells
  )
)


################################################################################
# STEP 1: Load
################################################################################

if (!file.exists(IN_PATH)) stop("Input not found: ", IN_PATH)

message("Loading RenalUrothelium_renal_epithelial.rds ...")
so <- readRDS(IN_PATH)
message(sprintf("  Loaded: %d cells × %d genes", ncol(so), nrow(so)))
log_mem("after load")

for (col in c("epithelial_compartment", "assign_method",
              "sample_id", "condition", "technology")) {
  if (col %in% colnames(so@meta.data)) {
    message(sprintf("  Cells per %s:", col))
    print(sort(table(so@meta.data[[col]]), decreasing = TRUE))
  }
}


################################################################################
# STEP 2: Set data layer = counts (already log-norm)
################################################################################

message("Setting data layer from counts (already log-norm) ...")
so[["RNA"]]$data <- so[["RNA"]]$counts
log_mem("after setting data layer")


################################################################################
# STEP 3: pct_mt + FindVariableFeatures
################################################################################

message("Computing pct_mt ...")
so[["pct_mt"]] <- PercentageFeatureSet(so, pattern = "^mt-")
message(sprintf("  pct_mt range: %.2f – %.2f", min(so$pct_mt), max(so$pct_mt)))

n_na <- sum(is.na(so$pct_mt))
if (n_na > 0) {
  med_val <- median(so$pct_mt, na.rm = TRUE)
  so$pct_mt[is.na(so$pct_mt)] <- med_val
  message(sprintf("  Imputed %d NAs in pct_mt with median (%.4f)", n_na, med_val))
}

message(sprintf("FindVariableFeatures (nfeatures = %d) ...", N_HVG))
so <- FindVariableFeatures(so, selection.method = "vst",
                           nfeatures = N_HVG, verbose = FALSE)
message(sprintf("  HVGs: %d", length(VariableFeatures(so))))


################################################################################
# STEP 4: ScaleData
################################################################################

message("ScaleData (regress pct_mt, HVGs only) ...")
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

so[["RNA"]]$counts     <- NULL
so[["RNA"]]$data       <- NULL
so[["RNA"]]$scale.data <- NULL
gc()
log_mem("after freeing expression data")


################################################################################
# STEP 6: RunHarmony
################################################################################

harmony_vars <- "sample_id"
for (v in c("technology", "source")) {
  if (v %in% colnames(so@meta.data) &&
      length(unique(so@meta.data[[v]])) > 1) {
    harmony_vars <- c(harmony_vars, v)
  }
}
message(sprintf("RunHarmony (batch = %s) ...",
                paste(harmony_vars, collapse = " + ")))

so <- RunHarmony(
  so,
  group.by.vars    = harmony_vars,
  reduction        = "pca",
  reduction.save   = "harmony",
  plot_convergence = FALSE,
  verbose          = FALSE
)
so[["pca"]] <- NULL
gc()
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

message("Generating UMAP overview plots ...")

make_dimplot <- function(grp, title, label = FALSE) {
  if (!grp %in% colnames(so@meta.data)) return(NULL)
  DimPlot(so, group.by = grp, reduction = "umap_harmony",
          label = label, repel = label, raster = TRUE) + ggtitle(title)
}

overview_plots <- Filter(Negate(is.null), list(
  make_dimplot("seurat_clusters",        "Harmony clusters",       label = TRUE),
  make_dimplot("epithelial_compartment", "Epithelial compartment", label = TRUE),
  make_dimplot("condition",              "By condition"),
  make_dimplot("sample_id",              "By sample"),
  make_dimplot("technology",             "By technology"),
  make_dimplot("assign_method",          "Assignment method"),
  make_dimplot("cell_type_original",     "Original cell type",     label = TRUE)
))

n_cols <- min(2L, length(overview_plots))
pdf(file.path(PLOT_DIR, "RenalEpithelial_UMAP_overview.pdf"),
    width = 18, height = ceiling(length(overview_plots) / n_cols) * 7)
print(wrap_plots(overview_plots, ncol = n_cols))
dev.off()
message("  Saved: RenalEpithelial_UMAP_overview.pdf")

# Per-compartment marker FeaturePlots — restore data layer temporarily
message("Restoring data layer for marker feature plots ...")
so_tmp <- readRDS(IN_PATH)
so[["RNA"]]$data <- so_tmp[["RNA"]]$counts
rm(so_tmp)
gc()

for (comp in names(MARKERS)) {
  genes   <- intersect(MARKERS[[comp]], rownames(so))
  missing <- setdiff(MARKERS[[comp]], rownames(so))
  message(sprintf("  %s: %d / %d markers present  (absent: %s)",
                  comp, length(genes), length(MARKERS[[comp]]),
                  if (length(missing)) paste(missing, collapse = ", ") else "none"))
  if (length(genes) == 0) next

  ncols <- min(5L, length(genes))
  nrows <- ceiling(length(genes) / ncols)
  pdf(file.path(PLOT_DIR, sprintf("RenalEpithelial_markers_%s.pdf", comp)),
      width = ncols * 4, height = nrows * 4)
  print(FeaturePlot(so, features = genes, reduction = "umap_harmony",
                    ncol = ncols, raster = TRUE))
  dev.off()
  message(sprintf("  Saved: RenalEpithelial_markers_%s.pdf", comp))
}

so[["RNA"]]$data <- NULL
gc()


################################################################################
# STEP 10: Save
################################################################################

message(sprintf("Saving to %s ...", OUT_PATH))
saveRDS(so, OUT_PATH)

message("\n===== Renal epithelial Harmony reclustering complete =====")
message(sprintf("Cells    : %s", format(ncol(so), big.mark = ",")))
message(sprintf("Clusters : %d", length(unique(so$seurat_clusters))))
message(sprintf("Output   : %s", OUT_PATH))
