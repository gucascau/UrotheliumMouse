################################################################################
# Script 02: QC, Filtering, and Doublet Removal
#
# Steps per sample:
#   1. Compute nFeature / nCount / % mitochondrial
#   2. Apply QC thresholds
#   3. Run DoubletFinder to predict and remove doublets
#   4. Save QC-filtered, doublet-free Seurat object
#
# Run AFTER 01_load_datasets.R has saved each _seurat.rds.
# Memory: reads one object at a time — 64 GB RAM is sufficient.
################################################################################

library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)
library(DoubletFinder)   # devtools::install_github("chris-mcginnis-ucsf/DoubletFinder")
library(patchwork)
library(future)          # parallelizes DoubletFinder paramSweep

# ── Parallelism setup ─────────────────────────────────────────────────────────
# doubletFinder(sct=TRUE) calls SCTransform internally on the artificial-doublet
# dataset. SCTransform and paramSweep both use future internally. future checks
# globals size even under plan(sequential), so the maxSizeOfObjects limit fires
# regardless of plan. Set Inf to disable the check; actual parallelism is
# controlled by BLAS/OMP threads (OMP_NUM_THREADS in the SLURM script).
plan(sequential)
options(future.globals.maxSize = Inf)
options(future.rng.onMisuse   = "ignore")
message("future plan: sequential | globals size limit: disabled")

# ── SLURM array index (0-based) ───────────────────────────────────────────────
# When submitted as an array job, each task processes one sample.
# When run sequentially (no --array-index flag), all samples are processed.
args        <- commandArgs(trailingOnly = TRUE)
array_idx   <- NULL
if ("--array-index" %in% args) {
  array_idx <- as.integer(args[which(args == "--array-index") + 1])
  message(sprintf("SLURM array mode: processing sample index %d", array_idx))
}

DATA_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
OBJ_DIR  <- file.path(DATA_DIR, "seurat_objects")
QC_DIR   <- file.path(DATA_DIR, "qc_plots")
dir.create(QC_DIR, showWarnings = FALSE)

################################################################################
# QC thresholds – adjust per dataset if needed
################################################################################

DEFAULT_QC <- list(
  min_features = 200,
  max_features = 8000,   # upper cutoff for first-pass (DoubletFinder refines this)
  max_mt_pct   = 40,     # % mitochondrial reads
  min_counts   = 300
)

# Dataset-specific overrides (name: list of thresholds)
# Set min_features = NULL to skip ALL QC+doublet removal for pre-processed datasets
QC_OVERRIDES <- list(
  EmbryosE9_5ToE13_5 = list(min_features = 100, max_features = 4000, max_mt_pct = 20,
                              min_counts = 200),   # sci-RNA-seq: lower UMI/cell
  KudoUUOUrothelium  = list(min_features = NULL),  # already QC'd + integrated
  MKA                = list(min_features = NULL),  # atlas reference – skip
  ChenSpatial        = list(min_features = NULL),  # spatial reference – skip
  LakesnRNA          = list(min_features = NULL),  # reference – skip
  KidneyUUO7         = list(min_features = NULL),  # pre-normalized floats: no counts layer, QC metrics undefined
  KidneyUUO8         = list(min_features = NULL)   # pre-normalized floats: no counts layer, QC metrics undefined
)

get_qc <- function(sid) {
  if (sid %in% names(QC_OVERRIDES)) {
    override <- QC_OVERRIDES[[sid]]
    base     <- DEFAULT_QC
    base[names(override)] <- override
    return(base)
  }
  DEFAULT_QC
}

################################################################################
# DoubletFinder parameters
################################################################################

# Hard cell-count limit: skip DoubletFinder for objects larger than this after
# QC filtering.  DoubletFinder is impractical above ~200k cells (memory + time).
DOUBLETFINDER_CELL_LIMIT <- 200000

# Datasets that already have doublet annotations in their source data and were
# filtered in 01_load_datasets.R — skip DoubletFinder even if below the limit.
SKIP_DOUBLETFINDER <- c(
  "EmbryosE9_5ToE13_5"   # sci-RNA-seq GSE119945: doublets annotated in cell_annotate.csv
)

# Expected doublet rate per dataset (10X Genomics ~ 0.8% per 1000 cells loaded)
# Adjust based on the library loading density of each sample.
# A simple heuristic: rate = 0.008 * (ncells / 1000), capped at ~0.25
expected_doublet_rate <- function(n_cells) {
  min(0.008 * (n_cells / 1000), 0.25)
}

# pN: proportion of artificial doublets (default 0.25)
# pK: neighborhood size – estimated per sample via paramSweep (see below)
DF_PN <- 0.25

################################################################################
# Helper: run DoubletFinder on a pre-processed Seurat object
################################################################################

run_doubletfinder <- function(so, doublet_rate, pN = DF_PN, dims = 1:15, sct_used = TRUE) {

  n_cells <- ncol(so)
  message(sprintf("    DoubletFinder: %d cells, expected doublet rate = %.1f%%",
                  n_cells, doublet_rate * 100))

  # --- 1. Estimate pK via parameter sweep ---
  sweep_res  <- paramSweep(so, PCs = dims, sct = sct_used, num.cores = 1)
  sweep_stat <- summarizeSweep(sweep_res, GT = FALSE)
  bcmvn      <- find.pK(sweep_stat)
  best_pK    <- as.numeric(as.character(
    bcmvn$pK[which.max(bcmvn$BCmetric)]
  ))
  message(sprintf("    Optimal pK = %.4f", best_pK))

  # --- 2. Estimate homotypic doublet proportion ---
  # Use seurat_clusters if available, otherwise skip homotypic correction
  if ("seurat_clusters" %in% colnames(so@meta.data)) {
    homotypic_prop <- modelHomotypic(so$seurat_clusters)
  } else {
    homotypic_prop <- 0
  }
  n_exp_poi       <- round(doublet_rate * n_cells)
  n_exp_poi_adj   <- round(n_exp_poi * (1 - homotypic_prop))
  message(sprintf("    nExp = %d  (adj = %d, homotypic prop = %.3f)",
                  n_exp_poi, n_exp_poi_adj, homotypic_prop))

  # --- 3. Run DoubletFinder ---
  so <- doubletFinder(so,
                      PCs        = dims,
                      pN         = pN,
                      pK         = best_pK,
                      nExp       = n_exp_poi_adj,
                      reuse.pANN = FALSE,
                      sct        = sct_used)

  # The classification column name is dynamic — find it
  df_col <- grep("^DF\\.classifications", colnames(so@meta.data), value = TRUE)
  df_col <- df_col[length(df_col)]   # take the last one
  so$doublet_class <- so@meta.data[[df_col]]

  pann_col <- grep("^pANN", colnames(so@meta.data), value = TRUE)
  pann_col <- pann_col[length(pann_col)]
  so$pANN <- so@meta.data[[pann_col]]

  return(so)
}


################################################################################
# Main loop: process each sample
################################################################################

# sort(method="radix") uses byte-order (locale-independent), matching bash LC_COLLATE=C.
# This ensures SLURM array task IDs map to the same samples regardless of locale.
rds_files  <- sort(
                list.files(OBJ_DIR, pattern = "_seurat\\.rds$", full.names = TRUE),
                method = "radix"
              )
sample_ids <- sub("_seurat\\.rds$", "", basename(rds_files))

# If running as a SLURM array job, process only the assigned sample
if (!is.null(array_idx)) {
  idx        <- array_idx + 1L   # R is 1-based
  rds_files  <- rds_files[idx]
  sample_ids <- sample_ids[idx]
  message(sprintf("Array mode: processing %s", sample_ids))
}

for (i in seq_along(rds_files)) {
  sid <- sample_ids[i]
  message(sprintf("\n=== QC + Doublet removal: %s ===", sid))

  so  <- readRDS(rds_files[i])
  qc  <- get_qc(sid)

  # ── Skip pre-processed/reference datasets ──────────────────────────────────
  if (is.null(qc$min_features)) {
    message("  Pre-processed dataset – skipping QC and DoubletFinder")
    saveRDS(so, file.path(OBJ_DIR, paste0(sid, "_qc.rds")))
    next
  }

  # ── Step 1: Mitochondrial percentage ──────────────────────────────────────
  so[["pct_mt"]] <- PercentageFeatureSet(so, pattern = "^mt-|^MT-")

  n_before <- ncol(so)
  message(sprintf("  Before QC: %d cells", n_before))

  # Violin plot before filtering
  p_before <- VlnPlot(so,
                       features = c("nFeature_RNA", "nCount_RNA", "pct_mt"),
                       ncol = 3, pt.size = 0) &
              theme(axis.title.x = element_blank())
  ggsave(file.path(QC_DIR, paste0(sid, "_qc_before.pdf")), p_before,
         width = 12, height = 4)

  # ── Step 2: Apply QC thresholds ───────────────────────────────────────────
  keep <- so$nFeature_RNA >= qc$min_features &
          so$nFeature_RNA <= qc$max_features &
          so$nCount_RNA   >= qc$min_counts   &
          so$pct_mt       <= qc$max_mt_pct

  so   <- so[, keep]
  message(sprintf("  After QC thresholds: %d cells (removed %d)",
                  ncol(so), n_before - ncol(so)))

  # ── Steps 3–5: Preprocessing + DoubletFinder (skipped for large / pre-annotated) ──
  skip_df <- sid %in% SKIP_DOUBLETFINDER || ncol(so) > DOUBLETFINDER_CELL_LIMIT
  if (skip_df) {
    if (sid %in% SKIP_DOUBLETFINDER) {
      message("  Skipping DoubletFinder: doublets already annotated/removed in source data")
    } else {
      message(sprintf("  Skipping DoubletFinder: %d cells exceeds limit of %d",
                      ncol(so), DOUBLETFINDER_CELL_LIMIT))
    }
  } else {
    prenorm <- "skip_normalization" %in% colnames(so@meta.data) && isTRUE(so$skip_normalization[1])

    if (prenorm) {
      # Pre-normalized data (GSE264184): values are already in the data layer.
      # Skip SCTransform; use FindVariableFeatures + ScaleData instead.
      message("  Pre-normalized data: skipping SCTransform (FindVariableFeatures + ScaleData)...")
      so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000,
                                 verbose = FALSE)
      so <- ScaleData(so, verbose = FALSE)
    } else {
      message("  Preprocessing for DoubletFinder (SCTransform + PCA + UMAP)...")
      # SCTransform uses future internally; switch to sequential to avoid
      # future.globals.maxSize errors. glmGamPoi (vst.flavor="v2") uses BLAS
      # threading and doesn't need future workers.
      plan(sequential)
      so <- SCTransform(so,
                        vars.to.regress = "pct_mt",
                        vst.flavor      = "v2",
                        verbose         = FALSE)
    }

    so <- RunPCA(so, npcs = 30, verbose = FALSE)
    so <- RunUMAP(so, dims = 1:15, verbose = FALSE)
    so <- FindNeighbors(so, dims = 1:15, verbose = FALSE)
    so <- FindClusters(so, resolution = 0.5, verbose = FALSE)

    dr <- expected_doublet_rate(ncol(so))
    so <- run_doubletfinder(so, doublet_rate = dr, dims = 1:15, sct_used = !prenorm)

    n_doublets <- sum(so$doublet_class == "Doublet")
    message(sprintf("  Doublets detected: %d (%.1f%%)",
                    n_doublets, 100 * n_doublets / ncol(so)))

    # Plot: UMAP colored by doublet classification
    p_umap <- DimPlot(so, group.by = "doublet_class",
                      cols = c(Singlet = "grey80", Doublet = "red")) +
              ggtitle(paste(sid, "– DoubletFinder")) +
              FeaturePlot(so, features = "pANN", reduction = "umap")
    ggsave(file.path(QC_DIR, paste0(sid, "_doublets_umap.pdf")), p_umap,
           width = 14, height = 6)

    so <- so[, so$doublet_class == "Singlet"]
    message(sprintf("  Final singlets: %d", ncol(so)))
  }

  # ── Step 6: Violin plot after all filtering ────────────────────────────────
  p_after <- VlnPlot(so,
                      features = c("nFeature_RNA", "nCount_RNA", "pct_mt"),
                      ncol = 3, pt.size = 0) &
             theme(axis.title.x = element_blank())
  ggsave(file.path(QC_DIR, paste0(sid, "_qc_after.pdf")), p_after,
         width = 12, height = 4)

  # ── Save ──────────────────────────────────────────────────────────────────
  saveRDS(so, file.path(OBJ_DIR, paste0(sid, "_qc.rds")))
  message(sprintf("  Saved: %s_qc.rds", sid))

  # Free memory before next sample
  rm(so); gc()
}

################################################################################
# Summary table
# Only written in sequential mode (no --array-index).
# In array mode, run this separately after all jobs finish:
#   Rscript 02_qc_and_filter.R --summary-only
################################################################################

if (is.null(array_idx) || "--summary-only" %in% args) {
  qc_files    <- list.files(OBJ_DIR, pattern = "_qc\\.rds$", full.names = TRUE)
  qc_ids      <- sub("_qc\\.rds$", "", basename(qc_files))
  cell_counts <- sapply(qc_files, function(f) ncol(readRDS(f)))

  summary_df  <- data.frame(sample_id = qc_ids, n_cells = cell_counts,
                             row.names = NULL)
  write.csv(summary_df, file.path(QC_DIR, "cell_counts_after_qc.csv"),
            row.names = FALSE)
  message("\nQC complete. Cell count summary:")
  print(summary_df)
}
message("Done: ", ifelse(is.null(array_idx), "all samples", sample_ids))
message("Filtered objects saved in: ", OBJ_DIR)
