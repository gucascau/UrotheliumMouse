################################################################################
# Script 04a: ScaleData + RunPCA  (pre-processing checkpoint)
#
# Reads the merged + log-normalized object from 03_integrate_harmony.R and runs:
#   1. Restore HVG selection
#   2. JoinLayers  — collapse per-sample layers into one matrix
#   3. ScaleData   — regress out pct_mt (HVGs only)
#   4. RunPCA      — 50 PCs on HVGs
#
# Writes:
#   integration_output/Integrated_ScaledPCA.rds
#
# Run BEFORE 04_integrate_harmony.R.
# Run AFTER  03_integrate_harmony.R has written:
#   integration_output/merged_normalized.rds
#   integration_output/hvg_list.rds
################################################################################

library(Seurat)

# Helper: log current R memory usage at key steps
log_mem <- function(label) {
  g  <- gc(verbose = FALSE)
  mb <- sum(g[, 2])
  message(sprintf("  [mem] %s: %.1f GB in use", label, mb / 1024))
}

# ── Paths ─────────────────────────────────────────────────────────────────────
DATA_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
OUT_DIR  <- file.path(DATA_DIR, "integration_output")
dir.create(OUT_DIR, showWarnings = FALSE)

NORM_PATH <- file.path(OUT_DIR, "merged_normalized.rds")
HVG_PATH  <- file.path(OUT_DIR, "hvg_list.rds")
PCA_PATH  <- file.path(OUT_DIR, "Integrated_ScaledPCA.rds")

# ── Parameters ────────────────────────────────────────────────────────────────
N_PCS <- 50


################################################################################
# STEP 1: Load merged normalized object + HVG list
################################################################################

for (p in c(NORM_PATH, HVG_PATH)) {
  if (!file.exists(p))
    stop("Required file not found: ", p,
         "\nRun 03_integrate_harmony.R first.")
}

message("Loading merged normalized object...")
merged <- readRDS(NORM_PATH)
message(sprintf("  Loaded: %d cells x %d genes", ncol(merged), nrow(merged)))
message(sprintf("  Layers: %s", paste(Layers(merged), collapse = ", ")))
log_mem("after load")

message("Loading HVG list...")
hvg_obj <- readRDS(HVG_PATH)
VariableFeatures(merged) <- hvg_obj$hvg
message(sprintf("  HVGs restored: %d", length(VariableFeatures(merged))))
rm(hvg_obj); gc()


################################################################################
# STEP 2: JoinLayers
################################################################################

message("JoinLayers — collapsing per-sample layers...")
merged <- JoinLayers(merged)
message(sprintf("  Layers after join: %s",
                paste(Layers(merged), collapse = ", ")))
log_mem("after JoinLayers")


################################################################################
# STEP 3: ScaleData
################################################################################

message("ScaleData (regress pct_mt, HVGs only)...")
merged <- ScaleData(merged,
                    features        = VariableFeatures(merged),
                    vars.to.regress = "pct_mt",
                    verbose         = FALSE)
message("  ScaleData done")
log_mem("after ScaleData")


################################################################################
# STEP 4: RunPCA
################################################################################

message(sprintf("RunPCA (%d PCs)...", N_PCS))
merged <- RunPCA(merged, npcs = N_PCS, verbose = FALSE)
message("  PCA done")
log_mem("after PCA")

message(sprintf("Saving PCA checkpoint to %s ...", PCA_PATH))
saveRDS(merged, PCA_PATH)

message("\n===== ScaleData + PCA complete =====")
message(sprintf("Total cells: %s", format(ncol(merged), big.mark = ",")))
message(sprintf("PCA dims:    %s", paste(dim(Embeddings(merged, "pca")),
                                         collapse = " x ")))
message(sprintf("Output:      %s", PCA_PATH))
