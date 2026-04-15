################################################################################
# Script 03c: Export merged_normalized.rds → AnnData (.h5ad)
#
# Reads the merged + log-normalised Seurat v5 object written by
# 03_integrate_harmony.R (merged_normalized.rds) and exports it as an h5ad
# file for downstream Python / scVI / scanpy workflows.
#
# AnnData layout:
#   adata.X              — log-normalised data (sparse float, cells × genes)
#   adata.layers["counts"] — raw integer counts, if present in merged object
#   adata.obs            — cell metadata (sample_id, condition, technology, …)
#   adata.var            — gene metadata + highly_variable flag (from hvg_list.rds)
#   adata.uns["hvg"]     — list of HVG gene names
#
# Memory: loading merged_normalized.rds + JoinLayers + transpose requires
#         ~400–512 GB RAM.  Run on a himem node (512 GB).
#
# Run AFTER 03_integrate_harmony.R has written:
#   integration_output/merged_normalized.rds
#   integration_output/hvg_list.rds
################################################################################

library(Seurat)
library(Matrix)
library(reticulate)

# ── Paths ─────────────────────────────────────────────────────────────────────
DATA_DIR  <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
OUT_DIR   <- file.path(DATA_DIR, "integration_output")
dir.create(OUT_DIR, showWarnings = FALSE)

NORM_PATH <- file.path(OUT_DIR, "merged_normalized.rds")
HVG_PATH  <- file.path(OUT_DIR, "hvg_list.rds")
H5AD_PATH <- file.path(OUT_DIR, "merged_normalized.h5ad")

message("=== 03c: Export merged_normalized.rds → h5ad ===")
message(sprintf("Input  : %s", NORM_PATH))
message(sprintf("Output : %s", H5AD_PATH))

for (p in c(NORM_PATH, HVG_PATH)) {
  if (!file.exists(p))
    stop(sprintf("Required file not found: %s\nRun 03_integrate_harmony.R first.", p))
}


################################################################################
# STEP 1: Load HVG list
################################################################################

message("\n[1/5] Loading HVG list...")
hvg_data  <- readRDS(HVG_PATH)
hvg_genes <- hvg_data$hvg
all_genes  <- hvg_data$all_genes
message(sprintf("  HVGs: %d / total genes: %d", length(hvg_genes), length(all_genes)))


################################################################################
# STEP 2: Load merged normalised Seurat object
################################################################################

message("\n[2/5] Loading merged_normalized.rds (this may take several minutes)...")
merged <- readRDS(NORM_PATH)
message(sprintf("  Loaded: %d cells × %d genes", ncol(merged), nrow(merged)))
message(sprintf("  Layers present: %s", paste(Layers(merged), collapse = ", ")))


################################################################################
# STEP 3: JoinLayers — collapse per-sample layers into single matrices
################################################################################

message("\n[3/5] JoinLayers — collapsing per-sample layers...")
DefaultAssay(merged) <- "RNA"
merged <- JoinLayers(merged)
message(sprintf("  Layers after join: %s", paste(Layers(merged), collapse = ", ")))


################################################################################
# STEP 4: Extract matrices and metadata
################################################################################

message("\n[4/5] Extracting data layer (log-normalised counts)...")

# Extract log-normalised matrix (genes × cells → will transpose below)
data_mat <- GetAssayData(merged, layer = "data")   # dgCMatrix: genes × cells
message(sprintf("  data layer: %d genes × %d cells  (class: %s)",
                nrow(data_mat), ncol(data_mat), class(data_mat)))

# Raw counts — present when all datasets had integer counts available.
# If counts layer is absent or differs in size, skip it to avoid errors.
counts_mat <- NULL
if ("counts" %in% Layers(merged)) {
  message("  Extracting counts layer (raw integer counts)...")
  counts_mat <- GetAssayData(merged, layer = "counts")
  if (!identical(dim(counts_mat), dim(data_mat))) {
    warning("  counts layer dimensions differ from data layer — skipping counts export")
    counts_mat <- NULL
  } else {
    message(sprintf("  counts layer: %d genes × %d cells", nrow(counts_mat), ncol(counts_mat)))
  }
} else {
  message("  No 'counts' layer found — exporting data layer only")
}

# Cell metadata — select standard columns, pad missing with NA
STANDARD_COLS <- c("sample_id", "condition", "technology", "paper", "gsm_id",
                   "species", "nFeature_RNA", "nCount_RNA", "pct_mt",
                   "cell_type_original", "cell_type", "doublet_class")
meta <- merged@meta.data[, intersect(colnames(merged@meta.data), STANDARD_COLS),
                         drop = FALSE]
for (col in setdiff(STANDARD_COLS, colnames(meta))) {
  meta[[col]] <- NA
}
meta <- meta[, STANDARD_COLS, drop = FALSE]
message(sprintf("  Cell metadata: %d cells × %d columns", nrow(meta), ncol(meta)))

# Gene (variable) metadata
var_df <- data.frame(
  gene_name       = rownames(data_mat),
  highly_variable = rownames(data_mat) %in% hvg_genes,
  row.names       = rownames(data_mat)
)
message(sprintf("  Gene metadata: %d genes (%d HVGs flagged)",
                nrow(var_df), sum(var_df$highly_variable)))

# Free the Seurat object before building Python objects
rm(merged); gc()


################################################################################
# STEP 5: Build AnnData and write h5ad
################################################################################

message("\n[5/5] Building AnnData and writing h5ad...")

anndata <- reticulate::import("anndata")
scipy   <- reticulate::import("scipy.sparse")

# Helper: convert a dgCMatrix (genes × cells) to scipy CSR (cells × genes)
dgc_to_csr <- function(mat, scipy_mod) {
  mat_t <- Matrix::t(mat)          # cells × genes (CSC)
  scipy_mod$csr_matrix(
    reticulate::tuple(
      mat_t@x,
      reticulate::tuple(mat_t@i, mat_t@p),
      reticulate::tuple(nrow(mat_t), ncol(mat_t))
    )
  )
}

message("  Converting data matrix to scipy CSR...")
X_csr <- dgc_to_csr(data_mat, scipy)
rm(data_mat); gc()

# Build AnnData with log-normalised data as X
adata <- anndata$AnnData(
  X   = X_csr,
  obs = meta,
  var = var_df
)

# Add raw counts as a layer if available
if (!is.null(counts_mat)) {
  message("  Converting counts matrix to scipy CSR...")
  counts_csr <- dgc_to_csr(counts_mat, scipy)
  rm(counts_mat); gc()
  adata$layers[["counts"]] <- counts_csr
  message("  Added adata.layers['counts']")
}

# Store HVG list in uns for downstream convenience
adata$uns[["hvg"]]   <- hvg_genes
adata$uns[["n_hvg"]] <- length(hvg_genes)

message(sprintf("\n  AnnData: %d cells × %d genes (%d HVGs)",
                as.integer(adata$n_obs), as.integer(adata$n_vars),
                length(hvg_genes)))
if ("sample_id" %in% colnames(meta)) {
  message(sprintf("  Samples: %s",
                  paste(sort(unique(meta$sample_id)), collapse = ", ")))
}

message(sprintf("\n  Writing to: %s", H5AD_PATH))
adata$write_h5ad(H5AD_PATH)

message("\n=== Export complete ===")
message(sprintf("Output: %s", H5AD_PATH))
