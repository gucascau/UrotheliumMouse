################################################################################
# Script 05: Export per-sample QC Seurat objects → individual h5ad files
#
# Reads every *_qc.rds in seurat_objects/ and writes one h5ad per sample to
# qc_h5ad/.  Samples already exported are skipped (restart-safe).
#
# AnnData layout per sample:
#   adata.X                — log-normalised data (if present), else raw counts
#   adata.layers["counts"] — raw integer counts (omitted when unavailable)
#   adata.obs              — cell metadata from Seurat @meta.data
#   adata.var              — gene names (+ highly_variable if hvg_list.rds exists)
#
# Ensembl ID → gene symbol conversion is applied automatically to datasets
# whose gene names start with ENSMUSG (ChenSpatial, LakesnRNA, MKA,
# KudoUUOUrothelium).
#
# Memory: individual samples are small; EmbryosE9_5ToE13_5 (~1.9 M cells) is
# the largest and needs ~128 GB.  Request 256 GB to be safe.
################################################################################

library(Seurat)
library(zellkonverter)
library(org.Mm.eg.db)

# ── Paths ─────────────────────────────────────────────────────────────────────
DATA_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
OBJ_DIR  <- file.path(DATA_DIR, "seurat_objects")
OUT_DIR  <- file.path(DATA_DIR, "qc_h5ad")
HVG_PATH <- file.path(DATA_DIR, "integration_output", "hvg_list.rds")

dir.create(OUT_DIR, showWarnings = FALSE)

message("=== 05: Export QC Seurat objects → h5ad ===")
message(sprintf("Input  : %s", OBJ_DIR))
message(sprintf("Output : %s", OUT_DIR))


################################################################################
# Helper: Ensembl ID → gene symbol (mouse, offline via org.Mm.eg.db)
# No-op when rownames are already gene symbols (none match ^ENSMUSG).
################################################################################

convert_ensembl_to_symbol <- function(so, sid) {
  genes <- rownames(so)
  if (!any(grepl("^ENSMUSG", genes))) return(so)

  n_ens <- sum(grepl("^ENSMUSG", genes))
  message(sprintf("    [Ensembl→Symbol] %s: %d Ensembl IDs detected (of %d genes)",
                  sid, n_ens, length(genes)))

  # Coerce v4 Assay → Assay5 if needed (ChenSpatial, LakesnRNA, MKA)
  if (!is(so[["RNA"]], "Assay5")) {
    message("    [Ensembl→Symbol] Coercing v4 Assay → Assay5")
    so[["RNA"]] <- as(so[["RNA"]], "Assay5")
  }
  so <- JoinLayers(so)

  sym_map <- setNames(
    suppressMessages(mapIds(org.Mm.eg.db,
                            keys      = genes,
                            column    = "SYMBOL",
                            keytype   = "ENSEMBL",
                            multiVals = "first")),
    genes
  )

  n_unmapped <- sum(is.na(sym_map))
  if (n_unmapped > 0) {
    message(sprintf("    [Ensembl→Symbol] %d genes unmapped — keeping original name",
                    n_unmapped))
    sym_map[is.na(sym_map)] <- names(sym_map)[is.na(sym_map)]
  }

  # Resolve duplicate symbols: keep highest-expressed gene per group
  dup_syms <- unique(sym_map[duplicated(sym_map)])
  if (length(dup_syms) > 0) {
    ref_layer  <- if ("counts" %in% Layers(so)) "counts" else "data"
    drop_genes <- character(0)
    for (sym in dup_syms) {
      grp  <- names(sym_map)[sym_map == sym]
      tots <- Matrix::rowSums(LayerData(so, layer = ref_layer)[grp, , drop = FALSE])
      drop_genes <- c(drop_genes, grp[order(tots, decreasing = TRUE)[-1]])
    }
    message(sprintf("    [Ensembl→Symbol] %d duplicate symbols resolved, %d genes dropped",
                    length(dup_syms), length(drop_genes)))
    so      <- so[!rownames(so) %in% drop_genes, ]
    sym_map <- sym_map[rownames(so)]
  }

  new_names <- unname(sym_map)
  for (lyr in Layers(so)) {
    mat           <- LayerData(so, layer = lyr)
    rownames(mat) <- new_names[match(rownames(mat), names(sym_map))]
    LayerData(so, layer = lyr) <- mat
  }
  rownames(so[["RNA"]]) <- new_names

  message(sprintf("    [Ensembl→Symbol] Done: %d genes retained", nrow(so)))
  so
}


################################################################################
# Load HVG list for highly_variable annotation (optional)
################################################################################

hvg_genes <- NULL
if (file.exists(HVG_PATH)) {
  hvg_genes <- readRDS(HVG_PATH)$hvg
  message(sprintf("HVG list loaded: %d genes", length(hvg_genes)))
} else {
  message("HVG list not found — skipping highly_variable annotation")
}


################################################################################
# Main loop: one h5ad per QC sample
################################################################################

qc_files   <- sort(list.files(OBJ_DIR, pattern = "_qc\\.rds$",
                               full.names = TRUE))
sample_ids <- sub("_qc\\.rds$", "", basename(qc_files))
n          <- length(qc_files)

message(sprintf("\nProcessing %d QC samples...\n", n))

for (i in seq_along(qc_files)) {
  sid      <- sample_ids[i]
  out_file <- file.path(OUT_DIR, paste0(sid, ".h5ad"))

  # ── Restart: skip already-exported samples ───────────────────────────────
  if (file.exists(out_file)) {
    message(sprintf("[%d/%d] %s — already exported, skipping", i, n, sid))
    next
  }

  message(sprintf("[%d/%d] %s", i, n, sid))

  # ── Load ─────────────────────────────────────────────────────────────────
  so <- readRDS(qc_files[i])

  # Ensure RNA is the active assay; rename if stored under a different name
  if (!"RNA" %in% Assays(so)) {
    cur <- DefaultAssay(so)
    message(sprintf("  Renaming assay '%s' → 'RNA'", cur))
    so  <- RenameAssays(so, assay.name = cur, new.assay.name = "RNA")
  }
  DefaultAssay(so) <- "RNA"

  # Drop non-RNA assays to keep the object lean
  for (assay in setdiff(Assays(so), "RNA")) so[[assay]] <- NULL

  # ── Ensembl → Symbol conversion ──────────────────────────────────────────
  so <- convert_ensembl_to_symbol(so, sid)

  # ── JoinLayers ───────────────────────────────────────────────────────────
  so <- JoinLayers(so)

  # ── Keep only data + counts (drop scale.data and any other extras) ───────
  keep_layers <- intersect(c("data", "counts"), Layers(so))
  so <- DietSeurat(so, layers = keep_layers, assays = "RNA")
  message(sprintf("  Layers exported: %s", paste(Layers(so), collapse = ", ")))
  message(sprintf("  Cells: %d  Genes: %d", ncol(so), nrow(so)))

  # ── Annotate var with highly_variable flag ───────────────────────────────
  if (!is.null(hvg_genes)) {
    so[["RNA"]]@meta.data[["highly_variable"]] <-
      rownames(so) %in% hvg_genes
  }

  # ── Convert to SCE and write h5ad ────────────────────────────────────────
  sce <- as.SingleCellExperiment(so)
  writeH5AD(sce, file = out_file)
  message(sprintf("  Saved: %s", out_file))

  rm(so, sce); gc()
}

message(sprintf("\n=== Export complete: h5ad files in %s ===", OUT_DIR))
