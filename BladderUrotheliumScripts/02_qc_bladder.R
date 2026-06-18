################################################################################
# Script 02: QC, Filtering, and Doublet Removal — Bladder Urothelium
#
# Steps per sample:
#   1. Compute nFeature / nCount / % mitochondrial
#   2. Apply QC thresholds
#   3. SCTransform + PCA + UMAP (preprocessing for DoubletFinder)
#   4. DoubletFinder — detect and remove doublets
#   5. Save QC-filtered, doublet-free Seurat object (_qc.rds)
#   6. Violin QC plots before/after filtering
#
# Supports SLURM array mode:
#   Rscript 02_qc_bladder.R --array-index <0-based-index>
# Sequential (all samples):
#   Rscript 02_qc_bladder.R
################################################################################

library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)
library(DoubletFinder)
library(patchwork)
library(future)

plan(sequential)
options(future.globals.maxSize = Inf)
options(future.rng.onMisuse   = "ignore")
message("future plan: sequential | globals size limit: disabled")

# ── SLURM array index ─────────────────────────────────────────────────────────
args      <- commandArgs(trailingOnly = TRUE)
array_idx <- NULL
if ("--array-index" %in% args) {
  array_idx <- as.integer(args[which(args == "--array-index") + 1])
  message(sprintf("SLURM array mode: processing sample index %d", array_idx))
}

DATA_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/BladderUrothelium"
OBJ_DIR  <- file.path(DATA_DIR, "seurat_objects")
QC_DIR   <- file.path(DATA_DIR, "qc_plots")
dir.create(QC_DIR, showWarnings = FALSE)

################################################################################
# QC thresholds
################################################################################

DEFAULT_QC <- list(
  min_features = 200,
  max_features = 8000,
  max_mt_pct   = 40,
  min_counts   = 300
)

# Dataset-specific overrides (set min_features = NULL to skip QC entirely)
QC_OVERRIDES <- list()

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

DOUBLETFINDER_CELL_LIMIT <- 200000
SKIP_DOUBLETFINDER       <- character(0)   # none pre-annotated in this dataset
DF_PN                    <- 0.25

expected_doublet_rate <- function(n_cells) min(0.008 * (n_cells / 1000), 0.25)

run_doubletfinder <- function(so, doublet_rate, pN = DF_PN,
                               dims = 1:15, sct_used = TRUE) {
  n_cells <- ncol(so)
  message(sprintf("    DoubletFinder: %d cells, expected rate = %.1f%%",
                  n_cells, doublet_rate * 100))

  sweep_res  <- paramSweep(so, PCs = dims, sct = sct_used, num.cores = 1)
  sweep_stat <- summarizeSweep(sweep_res, GT = FALSE)
  bcmvn      <- find.pK(sweep_stat)
  best_pK    <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  message(sprintf("    Optimal pK = %.4f", best_pK))

  homotypic_prop <- if ("seurat_clusters" %in% colnames(so@meta.data))
    modelHomotypic(so$seurat_clusters) else 0
  n_exp        <- round(doublet_rate * n_cells)
  n_exp_adj    <- round(n_exp * (1 - homotypic_prop))
  message(sprintf("    nExp = %d  (adj = %d, homotypic prop = %.3f)",
                  n_exp, n_exp_adj, homotypic_prop))

  so <- doubletFinder(so, PCs = dims, pN = pN, pK = best_pK,
                      nExp = n_exp_adj, reuse.pANN = FALSE, sct = sct_used)

  df_col   <- tail(grep("^DF\\.classifications", colnames(so@meta.data),
                        value = TRUE), 1)
  pann_col <- tail(grep("^pANN", colnames(so@meta.data), value = TRUE), 1)
  so$doublet_class <- so@meta.data[[df_col]]
  so$pANN          <- so@meta.data[[pann_col]]
  so
}


################################################################################
# Main processing loop
################################################################################

rds_files  <- sort(list.files(OBJ_DIR, pattern = "_seurat\\.rds$",
                               full.names = TRUE), method = "radix")
sample_ids <- sub("_seurat\\.rds$", "", basename(rds_files))

if (!is.null(array_idx)) {
  idx        <- array_idx + 1L
  rds_files  <- rds_files[idx]
  sample_ids <- sample_ids[idx]
  message(sprintf("Array mode: processing %s", sample_ids))
}

for (i in seq_along(rds_files)) {
  sid <- sample_ids[i]
  message(sprintf("\n=== QC + Doublet removal: %s ===", sid))

  so  <- readRDS(rds_files[i])
  qc  <- get_qc(sid)

  if (is.null(qc$min_features)) {
    message("  Pre-processed — skipping QC and DoubletFinder")
    saveRDS(so, file.path(OBJ_DIR, paste0(sid, "_qc.rds")))
    next
  }

  # ── Step 1: Mitochondrial percentage ──────────────────────────────────────
  so[["pct_mt"]] <- PercentageFeatureSet(so, pattern = "^mt-|^MT-")
  n_before <- ncol(so)
  message(sprintf("  Before QC: %d cells", n_before))

  p_before <- VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "pct_mt"),
                       ncol = 3, pt.size = 0) &
              theme(axis.title.x = element_blank())
  ggsave(file.path(QC_DIR, paste0(sid, "_qc_before.pdf")), p_before,
         width = 12, height = 4)

  # ── Step 2: Filter ────────────────────────────────────────────────────────
  keep <- so$nFeature_RNA >= qc$min_features &
          so$nFeature_RNA <= qc$max_features &
          so$nCount_RNA   >= qc$min_counts   &
          so$pct_mt       <= qc$max_mt_pct
  so   <- so[, keep]
  message(sprintf("  After QC thresholds: %d cells (removed %d)",
                  ncol(so), n_before - ncol(so)))

  # ── Steps 3–4: DoubletFinder (skip if too many cells or pre-annotated) ────
  skip_df <- sid %in% SKIP_DOUBLETFINDER || ncol(so) > DOUBLETFINDER_CELL_LIMIT
  if (skip_df) {
    if (sid %in% SKIP_DOUBLETFINDER) {
      message("  Skipping DoubletFinder: doublets already annotated in source")
    } else {
      message(sprintf("  Skipping DoubletFinder: %d cells exceeds limit of %d",
                      ncol(so), DOUBLETFINDER_CELL_LIMIT))
    }
  } else {
    message("  Preprocessing for DoubletFinder (SCTransform + PCA + UMAP)...")
    plan(sequential)
    so <- JoinLayers(so, assay = "RNA")   # Seurat v5: merge split layers before DoubletFinder
    so <- SCTransform(so, vars.to.regress = "pct_mt", vst.flavor = "v2",
                      verbose = FALSE)
    so <- RunPCA(so, npcs = 30, verbose = FALSE)
    so <- RunUMAP(so, dims = 1:15, verbose = FALSE)
    so <- FindNeighbors(so, dims = 1:15, verbose = FALSE)
    so <- FindClusters(so, resolution = 0.5, verbose = FALSE)

    dr <- expected_doublet_rate(ncol(so))
    so <- run_doubletfinder(so, doublet_rate = dr, dims = 1:15, sct_used = TRUE)

    n_doublets <- sum(so$doublet_class == "Doublet")
    message(sprintf("  Doublets detected: %d (%.1f%%)",
                    n_doublets, 100 * n_doublets / ncol(so)))

    p_umap <- DimPlot(so, group.by = "doublet_class",
                      cols = c(Singlet = "grey80", Doublet = "red")) +
              ggtitle(paste(sid, "– DoubletFinder")) +
              FeaturePlot(so, features = "pANN", reduction = "umap")
    ggsave(file.path(QC_DIR, paste0(sid, "_doublets_umap.pdf")), p_umap,
           width = 14, height = 6)

    so <- so[, so$doublet_class == "Singlet"]
    message(sprintf("  Final singlets: %d", ncol(so)))
  }

  # ── Step 5: Violin after filtering ────────────────────────────────────────
  p_after <- VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "pct_mt"),
                      ncol = 3, pt.size = 0) &
             theme(axis.title.x = element_blank())
  ggsave(file.path(QC_DIR, paste0(sid, "_qc_after.pdf")), p_after,
         width = 12, height = 4)

  saveRDS(so, file.path(OBJ_DIR, paste0(sid, "_qc.rds")))
  message(sprintf("  Saved: %s_qc.rds", sid))
  rm(so); gc()
}

################################################################################
# Summary (sequential mode or --summary-only)
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
