################################################################################
# Script 03b: Export QC-filtered count data to AnnData (.h5ad) for scVI/scANVI
#
# scVI requires RAW INTEGER COUNTS — log-normalised floats will bias the ELBO.
# Datasets that were deposited as pre-normalised matrices (no integer counts)
# cannot be used for scVI training and are excluded here:
#   - ChenSpatial    (SCE conversion, data layer only)
#   - LakesnRNA      (SCE conversion, data layer only)
#   - MKA            (SCE conversion, data layer only)
#   - KudoUUOUrothelium (already integrated RDS, data layer only)
#   - KidneyUUO7/8   (GSE264184 dense float matrix, data layer only)
#
# Output:
#   integration_output/scvi_input.h5ad
#     adata.X            — raw counts (sparse integer)
#     adata.obs          — cell metadata (sample_id, condition, technology, …)
#     adata.var          — gene metadata + highly_variable flag (from merged HVG)
#     adata.uns["hvg"]   — list of HVG gene names
#
# Run AFTER 03_integrate_harmony.R has written merged_normalized.rds
# (which carries the HVG selection).
################################################################################

library(Seurat)
library(Matrix)
library(reticulate)

# ── Paths ─────────────────────────────────────────────────────────────────────
DATA_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
OBJ_DIR  <- file.path(DATA_DIR, "seurat_objects")
OUT_DIR  <- file.path(DATA_DIR, "integration_output")
dir.create(OUT_DIR, showWarnings = FALSE)

H5AD_PATH <- file.path(OUT_DIR, "scvi_input.h5ad")

# ── Parameters ────────────────────────────────────────────────────────────────
N_HVG           <- 3000
EXCLUDE_SAMPLES <- c("BladderHomogenate1", "BladderHomogenate2")

# Samples that have ONLY a data layer (no raw counts).
# scVI requires integer counts; these are excluded from the h5ad.
NO_COUNTS_SAMPLES <- c(
  "ChenSpatial",         # SCE/zellkonvert conversion — only data layer
  "LakesnRNA",           # SCE/zellkonvert conversion — only data layer
  "MKA",                 # SCE/zellkonvert conversion — only data layer
  "KudoUUOUrothelium",   # pre-integrated RDS — data layer only
  "KidneyUUO7",          # GSE264184 dense float matrix
  "KidneyUUO8"           # GSE264184 dense float matrix
)

message("=== 03b: Export for scVI/scANVI ===")
message(sprintf("Excluded from integration : %s", paste(EXCLUDE_SAMPLES, collapse = ", ")))
message(sprintf("Excluded (no raw counts)  : %s", paste(NO_COUNTS_SAMPLES, collapse = ", ")))


################################################################################
# STEP 1: Load HVG list + full gene universe from merged_pre_integration.rds
################################################################################

hvg_path <- file.path(OUT_DIR, "hvg_list.rds")
if (!file.exists(hvg_path)) {
  stop(sprintf(
    "hvg_list.rds not found at: %s\nRun 03_integrate_harmony.R first.",
    hvg_path
  ))
}

message("Loading HVG list and gene universe from hvg_list.rds...")
hvg_data  <- readRDS(hvg_path)
hvg_genes <- hvg_data$hvg
all_genes <- hvg_data$all_genes
message(sprintf("  HVGs: %d / total genes: %d", length(hvg_genes), length(all_genes)))


################################################################################
# STEP 2: Load each QC object and collect counts + metadata
################################################################################

qc_files   <- sort(list.files(OBJ_DIR, pattern = "_qc\\.rds$", full.names = TRUE),
                   method = "radix")
sample_ids <- sub("_qc\\.rds$", "", basename(qc_files))

# Apply exclusions
keep <- !(sample_ids %in% EXCLUDE_SAMPLES | sample_ids %in% NO_COUNTS_SAMPLES)
qc_files   <- qc_files[keep]
sample_ids <- sample_ids[keep]

message(sprintf("\nProcessing %d count-based datasets:", length(sample_ids)))

counts_list <- vector("list", length(sample_ids))
meta_list   <- vector("list", length(sample_ids))

for (i in seq_along(qc_files)) {
  sid <- sample_ids[i]
  message(sprintf("  [%d/%d] %s", i, length(qc_files), sid))

  so <- readRDS(qc_files[i])
  DefaultAssay(so) <- "RNA"

  # Seurat v5: join per-sample layers into a single matrix before extraction
  so <- JoinLayers(so)

  # Confirm a counts layer exists
  if (!"counts" %in% Layers(so)) {
    warning(sprintf("  [SKIP] %s has no 'counts' layer after JoinLayers — skipping", sid))
    next
  }

  # Extract raw integer counts
  cnt <- GetAssayData(so, layer = "counts")   # dgCMatrix: genes × cells

  # Align gene universe to the shared gene list from the merged object.
  # Add zero rows for genes absent in this dataset.
  missing_genes <- setdiff(all_genes, rownames(cnt))
  if (length(missing_genes) > 0) {
    zero_rows <- Matrix::sparseMatrix(
      i = integer(0), j = integer(0),
      dims = c(length(missing_genes), ncol(cnt)),
      dimnames = list(missing_genes, colnames(cnt))
    )
    cnt <- rbind(cnt, zero_rows)
  }
  cnt <- cnt[all_genes, , drop = FALSE]   # reorder to shared gene universe

  # Make barcodes globally unique: prepend sample_id
  colnames(cnt) <- paste(sid, colnames(cnt), sep = "_")

  counts_list[[i]] <- cnt

  # Metadata — extract whichever standard columns are present, then pad the
  # rest with NA so every sample produces a data frame with identical columns.
  # rbind requires identical column sets across list entries.
  STANDARD_COLS <- c("sample_id", "condition", "technology", "paper", "gsm_id",
                     "species", "nFeature_RNA", "nCount_RNA", "pct_mt",
                     "cell_type_original", "cell_type", "doublet_class")
  meta <- so@meta.data[, intersect(colnames(so@meta.data), STANDARD_COLS),
                       drop = FALSE]
  for (col in setdiff(STANDARD_COLS, colnames(meta))) {
    meta[[col]] <- NA
  }
  meta <- meta[, STANDARD_COLS, drop = FALSE]   # enforce consistent column order
  rownames(meta) <- colnames(cnt)
  meta_list[[i]] <- meta

  rm(so, cnt); gc()
}

# Remove any skipped entries
counts_list <- Filter(Negate(is.null), counts_list)
meta_list   <- Filter(Negate(is.null), meta_list)


################################################################################
# STEP 3: Concatenate into one sparse matrix
################################################################################

message("\nConcatenating count matrices...")
counts_all <- do.call(cbind, counts_list)
meta_all   <- do.call(rbind, meta_list)
rm(counts_list, meta_list); gc()

message(sprintf("Combined: %d cells × %d genes", ncol(counts_all), nrow(counts_all)))

# Fill NA pct_mt with 0 for pre-processed datasets
if ("pct_mt" %in% colnames(meta_all)) {
  na_mt <- is.na(meta_all$pct_mt)
  if (any(na_mt)) {
    message(sprintf("Filling %d NA pct_mt values with 0", sum(na_mt)))
    meta_all$pct_mt[na_mt] <- 0
  }
}


################################################################################
# STEP 4: Build and write AnnData via reticulate
################################################################################

message("Writing h5ad via reticulate + anndata...")

# Ensure anndata is importable in the active Python environment.
# Set RETICULATE_PYTHON or use use_condaenv() before running if needed.
anndata <- reticulate::import("anndata")
scipy   <- reticulate::import("scipy.sparse")

# AnnData expects cells × genes (transpose of R's genes × cells)
counts_t <- Matrix::t(counts_all)   # cells × genes, still sparse

# Convert to scipy CSR for efficient row-slicing in Python
counts_csr <- scipy$csr_matrix(
  reticulate::tuple(
    counts_t@x,
    reticulate::tuple(counts_t@i, counts_t@p),
    reticulate::tuple(nrow(counts_t), ncol(counts_t))
  )
)

# Variable (gene) metadata
var_df <- data.frame(
  gene_name      = all_genes,
  highly_variable = all_genes %in% hvg_genes,
  row.names       = all_genes
)

# Build AnnData
adata <- anndata$AnnData(
  X   = counts_csr,
  obs = meta_all,
  var = var_df
)

# Store HVG list and excluded samples in uns for reference
adata$uns[["hvg"]]              <- hvg_genes
adata$uns[["excluded_samples"]] <- c(EXCLUDE_SAMPLES, NO_COUNTS_SAMPLES)
adata$uns[["n_hvg"]]            <- length(hvg_genes)

message(sprintf("AnnData: %d cells × %d genes (%d HVGs)",
                as.integer(adata$n_obs), as.integer(adata$n_vars),
                length(hvg_genes)))
message(sprintf("Samples included: %s", paste(sort(unique(meta_all$sample_id)), collapse = ", ")))

adata$write_h5ad(H5AD_PATH)
message(sprintf("\nSaved: %s", H5AD_PATH))
message("=== Export complete ===")
