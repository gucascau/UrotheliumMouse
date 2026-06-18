################################################################################
# Script 03: Export per-sample QC Seurat objects → individual h5ad files
#
# Reads every *_qc.rds in seurat_objects/ and writes one h5ad per sample to
# qc_h5ad/. Samples already exported are skipped (restart-safe).
#
# AnnData layout per sample:
#   adata.X                — log-normalised data (if present), else raw counts
#   adata.layers["counts"] — raw integer counts
#   adata.obs              — cell metadata from Seurat @meta.data
#   adata.var              — gene names
#
# All bladder samples use gene symbols (no Ensembl conversion needed).
# Gene names were already set to symbols during loading (01_load_bladder.R).
################################################################################

library(Seurat)
library(zellkonverter)

DATA_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/BladderUrothelium"
OBJ_DIR  <- file.path(DATA_DIR, "seurat_objects")
OUT_DIR  <- file.path(DATA_DIR, "qc_h5ad")
dir.create(OUT_DIR, showWarnings = FALSE)

message("=== 03: Export QC Seurat objects → h5ad ===")
message(sprintf("Input  : %s", OBJ_DIR))
message(sprintf("Output : %s", OUT_DIR))

qc_files   <- sort(list.files(OBJ_DIR, pattern = "_qc\\.rds$",
                               full.names = TRUE))
sample_ids <- sub("_qc\\.rds$", "", basename(qc_files))
n          <- length(qc_files)

message(sprintf("\nProcessing %d QC samples...\n", n))

for (i in seq_along(qc_files)) {
  sid      <- sample_ids[i]
  out_file <- file.path(OUT_DIR, paste0(sid, ".h5ad"))

  if (file.exists(out_file)) {
    message(sprintf("[%d/%d] %s — already exported, skipping", i, n, sid))
    next
  }

  message(sprintf("[%d/%d] %s", i, n, sid))
  so <- readRDS(qc_files[i])

  # Ensure RNA is the active assay
  if (!"RNA" %in% Assays(so)) {
    cur <- DefaultAssay(so)
    message(sprintf("  Renaming assay '%s' → 'RNA'", cur))
    so  <- RenameAssays(so, assay.name = cur, new.assay.name = "RNA")
  }
  DefaultAssay(so) <- "RNA"

  # Drop any non-RNA assays (e.g. SCT from DoubletFinder preprocessing)
  for (assay in setdiff(Assays(so), "RNA")) {
    so[[assay]] <- NULL
  }

  so <- JoinLayers(so)

  # Keep only counts and data layers
  keep_layers <- intersect(c("data", "counts"), Layers(so))
  so <- DietSeurat(so, layers = keep_layers, assays = "RNA")
  message(sprintf("  Layers: %s  |  Cells: %d  Genes: %d",
                  paste(Layers(so), collapse = ", "), ncol(so), nrow(so)))

  sce <- as.SingleCellExperiment(so)
  writeH5AD(sce, file = out_file)
  message(sprintf("  Saved: %s", out_file))

  rm(so, sce); gc()
}

message(sprintf("\n=== Export complete: %s ===", OUT_DIR))
