################################################################################
# 04a_export_for_scvi.R
#
# Exports AllUrothelium_markers_gated.rds into the sparse-matrix format
# expected by 04b_scvi_AllUrothelium.py:
#
#   output/AllUrothelium_scvi_input/
#     matrix.mtx   – genes × cells log-normalised sparse matrix (HVGs only)
#     barcodes.tsv – cell barcodes
#     features.tsv – gene names  (matching matrix rows)
#     metadata.csv – full cell metadata
#
# Why log-normalised (data layer) and not raw counts?
#   The scVI-source files (bladder, kidney in-vivo) only have scVI
#   log-normalised expression stored as their "counts" layer — the original
#   raw UMI counts are not available.  The qc-source files have true raw counts,
#   but scVI requires a single likelihood for all cells.  We therefore export
#   the log-normalised "data" layer for every cell and train scVI with
#   gene_likelihood = "normal" (Gaussian decoder), which is the documented
#   approach for pre-normalised input.
#
# HVG strategy:
#   HVGs are computed per-sample independently (batch-aware) then ranked by
#   cross-sample frequency, preventing chemistry / organoid-specific genes
#   from dominating the feature set before scVI even sees the data.
#
# Input : output/AllUrothelium_markers_gated.rds
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

log_mem <- function(label) {
  g  <- gc(verbose = FALSE)
  mb <- sum(g[, 2])
  message(sprintf("  [mem] %s: %.1f GB in use", label, mb / 1024))
}

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/FinalUrotheliumCells"
SCR_DIR  <- file.path(BASE_DIR, "UrotheliumScripts")
OUT_DIR  <- file.path(SCR_DIR, "output")
SCVI_IN  <- file.path(OUT_DIR, "AllUrothelium_scvi_input")

args     <- commandArgs(trailingOnly = TRUE)
RDS_PATH <- if (length(args) >= 1) args[[1]] else
  file.path(OUT_DIR, "AllUrothelium_markers_gated.rds")

dir.create(SCVI_IN, showWarnings = FALSE, recursive = TRUE)

# ── Parameters ────────────────────────────────────────────────────────────────
N_HVG     <- 4000L   # genes to export; scVI can use all of them
MIN_CELLS <- 20L     # min cells per sample to include in HVG computation

# ── Load ──────────────────────────────────────────────────────────────────────
message("Loading: ", RDS_PATH)
if (!file.exists(RDS_PATH)) stop("RDS not found: ", RDS_PATH)
so <- readRDS(RDS_PATH)
message(sprintf("  %s cells × %s genes",
                format(ncol(so), big.mark = ","),
                format(nrow(so), big.mark = ",")))
log_mem("after load")

# Ensure a single joined data layer is present
if (!"data" %in% Layers(so[["RNA"]])) {
  message("  Joining layers to restore data layer ...")
  so <- JoinLayers(so)
}

# ── Batch-aware HVG selection ─────────────────────────────────────────────────
# Compute HVGs per sample, then rank genes by how many samples called them
# variable.  Genes that are only highly variable in one chemistry / condition
# drop to the bottom.
message(sprintf("Computing batch-aware HVGs (target %d, min_cells_per_sample = %d) ...",
                N_HVG, MIN_CELLS))

sample_ids <- unique(so$sample_id)
hvg_freq   <- setNames(integer(nrow(so)), rownames(so))
n_valid    <- 0L

for (sid in sample_ids) {
  idx <- which(so$sample_id == sid)
  if (length(idx) < MIN_CELLS) {
    message(sprintf("  Skipping %s (%d cells < %d)", sid, length(idx), MIN_CELLS))
    next
  }
  tryCatch({
    so_s  <- so[, idx]
    so_s  <- FindVariableFeatures(so_s, selection.method = "vst",
                                  nfeatures = N_HVG, verbose = FALSE)
    top_g <- VariableFeatures(so_s)
    hvg_freq[top_g] <- hvg_freq[top_g] + 1L
    n_valid <- n_valid + 1L
    rm(so_s)
  }, error = function(e) {
    message(sprintf("  HVG failed for %s: %s", sid, conditionMessage(e)))
  })
}

hvg <- names(sort(hvg_freq, decreasing = TRUE))[seq_len(N_HVG)]
message(sprintf("  %d HVGs from %d valid samples.", length(hvg), n_valid))
message(sprintf("  Top 10: %s", paste(head(hvg, 10), collapse = ", ")))
log_mem("after HVG")

# ── Extract data layer for selected HVGs ──────────────────────────────────────
message("Extracting data layer ...")
data_mat <- GetAssayData(so, assay = "RNA", layer = "data")[hvg, , drop = FALSE]
# data_mat is genes × cells (standard R sparse orientation)
pct_nz <- 100 * nnzero(data_mat) / prod(dim(data_mat))
message(sprintf("  Matrix: %d genes × %d cells  (%.1f%% non-zero)",
                nrow(data_mat), ncol(data_mat), pct_nz))

# ── Write sparse matrix (genes × cells) ───────────────────────────────────────
message("Writing matrix.mtx ...")
writeMM(data_mat, file.path(SCVI_IN, "matrix.mtx"))

message("Writing barcodes.tsv ...")
writeLines(colnames(data_mat), file.path(SCVI_IN, "barcodes.tsv"))

message("Writing features.tsv ...")
writeLines(rownames(data_mat), file.path(SCVI_IN, "features.tsv"))

# ── Metadata ──────────────────────────────────────────────────────────────────
message("Writing metadata.csv ...")
meta <- so@meta.data

# Fill NA in columns used as scVI covariates (scVI errors on NA covariates)
for (col in c("sample_id", "technology", "LabelClass")) {
  if (col %in% colnames(meta) && any(is.na(meta[[col]])))
    meta[[col]][is.na(meta[[col]])] <- "Unknown"
}
if ("pct_mt" %in% colnames(meta) && any(is.na(meta$pct_mt)))
  meta$pct_mt[is.na(meta$pct_mt)] <- median(meta$pct_mt, na.rm = TRUE)

write.csv(meta, file.path(SCVI_IN, "metadata.csv"))

# ── Print summary ─────────────────────────────────────────────────────────────
message(sprintf("\nExport complete → %s", SCVI_IN))
message(sprintf("  Cells       : %s", format(ncol(so), big.mark = ",")))
message(sprintf("  HVGs        : %d", length(hvg)))
message(sprintf("  Samples     : %d", length(unique(so$sample_id))))
if ("technology" %in% colnames(meta))
  message(sprintf("  Technologies: %s",
                  paste(sort(unique(meta$technology[!is.na(meta$technology)])),
                        collapse = ", ")))
if ("LabelClass" %in% colnames(meta))
  message(sprintf("  Label Classes     : %s",
                  paste(sort(unique(meta$LabelClass)), collapse = ", ")))
if ("Categories" %in% colnames(meta))
  message(sprintf("  Categories  : %s",
                  paste(sort(unique(meta$Categories)), collapse = ", ")))
