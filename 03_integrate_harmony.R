################################################################################
# Script 03: Merge and Integrate with Harmony (each sample as a batch)
#
# Workflow:
#   1. Load all QC-filtered Seurat objects
#   2. Merge into one object (Seurat v5 keeps per-sample layers)
#   3. Normalize per layer (NormalizeData — recommended before Harmony)
#   4. FindVariableFeatures across all samples
#   5. ScaleData (regress pct_mt)
#   6. RunPCA
#   7. RunHarmony — group.by.vars = "sample_id" (each sample = one batch)
#   8. FindNeighbors + FindClusters + RunUMAP on Harmony embedding
#
# For >500k cells: sketch-based workflow (SKETCH_CELLS) subsamples each
# sample for PCA/Harmony, then projects back to all cells.
#
# Memory: request himem (512 GB) for full dataset; sketch needs ~256 GB.
################################################################################

library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)
library(patchwork)
library(org.Mm.eg.db)   # offline Ensembl → gene symbol mapping for mouse

DATA_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
OBJ_DIR  <- file.path(DATA_DIR, "seurat_objects")
OUT_DIR  <- file.path(DATA_DIR, "integration_output")
dir.create(OUT_DIR, showWarnings = FALSE)

################################################################################
# Helper: convert Ensembl IDs → gene symbols (mouse, offline via org.Mm.eg.db)
#
# Called per-sample inside the load loop. No-op when rownames are already
# gene symbols (i.e. none match ^ENSMUSG).
#
# Duplicate-symbol resolution: for each group of Ensembl IDs that map to the
# same symbol, keep the one with the highest total counts/expression and drop
# the rest (typically pseudogenes or low-quality isoforms).
# Unmapped IDs (no symbol in org.Mm.eg.db) keep their Ensembl ID as a fallback
# so no genes are silently lost.
################################################################################

convert_ensembl_to_symbol <- function(so, sid) {
  genes <- rownames(so)

  # ── 1. Detect: skip if no ENSMUSG IDs present ────────────────────────────────
  if (!any(grepl("^ENSMUSG", genes))) return(so)
  n_ens <- sum(grepl("^ENSMUSG", genes))
  message(sprintf("    [Ensembl→Symbol] %s: %d Ensembl IDs detected (of %d genes)",
                  sid, n_ens, length(genes)))

  # ── 2. Upgrade v4 Assay → Assay5 so we can use a single layer-based code path
  # ChenSpatial / LakesnRNA / MKA were saved as v4 Assay objects (zellkonvert).
  # UpdateSeuratObject only bumps version metadata — it does NOT change the class.
  # as(assay, "Assay5") is the correct coercion that migrates slots → layers.
  if (!is(so[["RNA"]], "Assay5")) {
    message("    [Ensembl→Symbol] Coercing v4 Assay → Assay5")
    so[["RNA"]] <- as(so[["RNA"]], "Assay5")
  }

  # Collapse any per-sample layers into one matrix for uniform access
  so <- JoinLayers(so)

  # ── 3. Build the symbol map ──────────────────────────────────────────────────
  # mapIds with keytype="ENSEMBL":
  #   - ENSMUSG* genes  → mapped symbol or NA
  #   - gene symbols already in the list → NA (not valid Ensembl IDs)
  # This handles the mixed case (KudoUUOUrothelium has both in its gene list).
  sym_map <- setNames(
    suppressMessages(mapIds(org.Mm.eg.db,
                            keys      = genes,
                            column    = "SYMBOL",
                            keytype   = "ENSEMBL",
                            multiVals = "first")),
    genes
  )

  # ── 4. Fallback for unmapped / already-symbol genes ──────────────────────────
  n_unmapped <- sum(is.na(sym_map))
  if (n_unmapped > 0) {
    message(sprintf("    [Ensembl→Symbol] %d genes unmapped — keeping original name",
                    n_unmapped))
    sym_map[is.na(sym_map)] <- names(sym_map)[is.na(sym_map)]
  }

  # ── 5. Resolve duplicate symbols ─────────────────────────────────────────────
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

  # ── 6. Rename features in every layer + update the Assay5 feature index ──────
  for (lyr in Layers(so)) {
    mat           <- LayerData(so, layer = lyr)
    rownames(mat) <- new_names[match(rownames(mat), names(sym_map))]
    LayerData(so, layer = lyr) <- mat
  }
  rownames(so[["RNA"]]) <- new_names

  message(sprintf("    [Ensembl→Symbol] Done: %d genes retained", nrow(so)))
  so
}

# ── Parameters ────────────────────────────────────────────────────────────────
N_HVG        <- 3000    # highly variable genes
N_PCS        <- 50      # PCs to compute
HARMONY_DIMS <- 1:50    # Harmony dims to use downstream
RESOLUTION   <- 0.5     # clustering resolution

# Harmony batch variables — correct for per-sample effects AND cross-technology
# platform differences (Drop-seq / 10X v2 / 10X v3 / sci-RNA-seq / snRNA-seq / Multiome)
HARMONY_VARS <- c("sample_id", "technology")

# Sketch workflow is disabled: SketchData's LeverageScore method fails on
# mixed-layer merged objects (pre-processed atlas layers have a subset gene
# universe that doesn't contain all HVGs).  Run on the full dataset instead.
# 512 GB RAM is sufficient: ScaleData on 3000 HVGs × 2.8 M cells ≈ 70 GB dense.
SKETCH_CELLS <- NULL

# Samples to exclude from integration (e.g. non-kidney tissue)
EXCLUDE_SAMPLES <- c("BladderHomogenate1", "BladderHomogenate2")   # bladder homogenate — not kidney
# ──────────────────────────────────────────────────────────────────────────────


################################################################################
# STEP 1: Load QC-filtered objects from seurat_objects/
################################################################################

qc_files   <- list.files(OBJ_DIR, pattern = "_qc\\.rds$", full.names = TRUE)
sample_ids <- sub("_qc\\.rds$", "", basename(qc_files))

keep_idx   <- !sample_ids %in% EXCLUDE_SAMPLES
qc_files   <- qc_files[keep_idx]
sample_ids <- sample_ids[keep_idx]

if (length(EXCLUDE_SAMPLES) > 0)
  message(sprintf("Excluding: %s", paste(EXCLUDE_SAMPLES, collapse = ", ")))
message(sprintf("Loading %d QC-filtered datasets...", length(qc_files)))

seurat_list <- lapply(seq_along(qc_files), function(i) {
  sid <- sample_ids[i]
  message(sprintf("  [%d/%d] %s", i, length(qc_files), sid))
  so  <- readRDS(qc_files[i])

  if (!"RNA" %in% Assays(so)) {
    cur <- DefaultAssay(so)
    message(sprintf("    Renaming assay '%s' → 'RNA' for %s", cur, sid))
    so  <- RenameAssays(so, assay.name = cur, new.assay.name = "RNA")
  }
  DefaultAssay(so) <- "RNA"

  for (assay in c("SCT", "ATAC", "Spatial", "sketch", "originalexp")) {
    if (assay %in% Assays(so)) so[[assay]] <- NULL
  }

  # Convert Ensembl IDs → gene symbols before merging so all samples share
  # the same gene namespace. No-op for samples already using gene symbols.
  # Affected: ChenSpatial, LakesnRNA, MKA (all ENSMUSG), KudoUUOUrothelium (mixed).
  so <- convert_ensembl_to_symbol(so, sid)

  so <- RenameCells(so, add.cell.id = sid)
  so
})
names(seurat_list) <- sample_ids
message(sprintf("Total cells loaded: %s",
                format(sum(sapply(seurat_list, ncol)), big.mark = ",")))


################################################################################
# STEP 2: Merge — Seurat v5 keeps one layer per sample
################################################################################

message("Merging all objects (Seurat v5 layer-aware merge)...")
merged <- merge(
  x          = seurat_list[[1]],
  y          = seurat_list[-1],
  add.cell.ids = NULL,
  merge.data   = TRUE
)
rm(seurat_list); gc()

message(sprintf("Merged object: %s cells, %s genes",
                format(ncol(merged), big.mark = ","),
                format(nrow(merged), big.mark = ",")))
message(sprintf("Layers: %s", paste(Layers(merged), collapse = ", ")))

# ── Ensure pct_mt for all cells ──────────────────────────────────────────────
if (!"pct_mt" %in% colnames(merged@meta.data)) {
  message("pct_mt not found — computing for all cells...")
  merged[["pct_mt"]] <- PercentageFeatureSet(merged, pattern = "^mt-|^MT-")
} else {
  n_missing <- sum(is.na(merged$pct_mt))
  if (n_missing > 0) {
    message(sprintf("Filling %s missing pct_mt values with 0 (pre-processed)",
                    format(n_missing, big.mark = ",")))
    merged$pct_mt[is.na(merged$pct_mt)] <- 0
  }
}

batch_tbl <- table(merged$sample_id)
message("\nCells per sample (batch):")
for (nm in names(batch_tbl))
  message(sprintf("  %-20s : %s", nm, format(batch_tbl[nm], big.mark = ",")))


################################################################################
# STEP 3: Normalise per layer
################################################################################

message("\nNormalizing (log-normalization per sample layer)...")
merged <- NormalizeData(merged, normalization.method = "LogNormalize",
                        scale.factor = 10000, verbose = FALSE)


################################################################################
# STEP 4: Find highly variable genes
################################################################################

message(sprintf("Finding %d variable features...", N_HVG))
merged <- FindVariableFeatures(merged, nfeatures = N_HVG,
                               selection.method = "vst", verbose = FALSE)
message(sprintf("  Top HVGs: %s",
                paste(head(VariableFeatures(merged), 10), collapse = ", ")))

# Save tiny sidecar for 03b_export_for_scvi.R (HVG list + gene universe)
saveRDS(list(hvg = VariableFeatures(merged), all_genes = rownames(merged)),
        file.path(OUT_DIR, "hvg_list.rds"))
message("  HVG list saved (hvg_list.rds)")

# Save merged normalized object for downstream integration step (04_integrate_harmony.R)
norm_path <- file.path(OUT_DIR, "merged_normalized.rds")
message(sprintf("Saving merged normalized object to %s ...", norm_path))
saveRDS(merged, norm_path)
message("  merged_normalized.rds saved")

message("\n===== Prep complete (STEPS 1-4) =====")
message(sprintf("Total cells : %s", format(ncol(merged), big.mark = ",")))
message(sprintf("HVGs saved  : %d  →  %s", length(VariableFeatures(merged)),
                file.path(OUT_DIR, "hvg_list.rds")))
message(sprintf("Merged obj  : %s", norm_path))
message("Next: run 04_integrate_harmony.R for ScaleData / PCA / Harmony / UMAP")

################################################################################
# STEP 5: Sketch-based or full integration
################################################################################


if (!is.null(SKETCH_CELLS)) {

  ## ── Sketch workflow ─────────────────────────────────────────────────────────
  message(sprintf("\nSketch-based integration: %d cells per sample...", SKETCH_CELLS))

  # LeverageScore iterates per layer and fails when pre-processed atlas layers
  # (ChenSpatial, LakesnRNA, MKA — from SCE conversion) have a subset gene
  # universe that doesn't include all HVGs.  Uniform sampling is robust to
  # mixed-layer objects and is sufficient for a 2.8 M-cell dataset.
  merged <- SketchData(
    object         = merged,
    ncells         = SKETCH_CELLS,
    method         = "LeverageScore",
    sketched.assay = "sketch"
  )
  DefaultAssay(merged) <- "sketch"

  message(sprintf("Sketch size: %d cells", ncol(merged[["sketch"]])))

  # The sketch assay inherits per-sample layers from RNA.  ScaleData calls
  # StitchMatrix -> LayerData per layer and hits the same "features not found"
  # error as LeverageScore did.  JoinLayers collapses them into one matrix so
  # ScaleData sees a single layer with the full gene universe.
  merged <- JoinLayers(merged, assay = "sketch")

  # Scale and PCA on sketch
  message("ScaleData + RunPCA on sketch...")

  merged <- ScaleData(merged,
                      vars.to.regress = "pct_mt",
                      verbose = FALSE)
  merged <- RunPCA(merged, npcs = N_PCS, verbose = FALSE)

  # Harmony on sketch — each sample_id is one batch
  message(sprintf("RunHarmony on sketch (batch = %s)...", HARMONY_VARS))
  merged <- RunHarmony(
    object         = merged,
    group.by.vars  = HARMONY_VARS,
    reduction      = "pca",
    reduction.save = "harmony",
    max_iter       = 20,
    verbose        = FALSE
  )

  merged <- FindNeighbors(merged, reduction = "harmony", dims = HARMONY_DIMS,
                          verbose = FALSE)
  merged <- FindClusters(merged, resolution = RESOLUTION, verbose = FALSE)
  merged <- RunUMAP(merged, reduction = "harmony", dims = HARMONY_DIMS,
                    reduction.name = "umap.sketch", verbose = FALSE)

  # Project Harmony embedding to all cells
  message("Projecting Harmony embedding to full dataset...")
  merged <- ProjectIntegration(
    object         = merged,
    sketched.assay = "sketch",
    assay          = "RNA",
    reduction      = "harmony"
  )
  merged <- ProjectData(
    object             = merged,
    sketched.assay     = "sketch",
    assay              = "RNA",
    sketched.reduction = "harmony",
    full.reduction     = "full.harmony",
    dims               = HARMONY_DIMS,
    refdata            = list(seurat_clusters = "seurat_clusters")
  )
  merged <- RunUMAP(merged, reduction = "full.harmony", dims = HARMONY_DIMS,
                    reduction.name = "umap", verbose = FALSE)

} else {

  ## ── Full-dataset workflow ───────────────────────────────────────────────────
  message("\nFull-dataset integration (no sketch)...")

  # JoinLayers collapses per-sample layers into one matrix for ScaleData/PCA
  merged <- JoinLayers(merged)

  message("ScaleData + RunPCA on full dataset...")
  merged <- ScaleData(merged,
                      vars.to.regress = "pct_mt",
                      verbose = FALSE)
  merged <- RunPCA(merged, npcs = N_PCS, verbose = FALSE)

  # Harmony — each sample_id is one batch
  message(sprintf("RunHarmony on full dataset (batch = %s)...", HARMONY_VARS))
  merged <- RunHarmony(
    object         = merged,
    group.by.vars  = HARMONY_VARS,
    reduction      = "pca",
    reduction.save = "harmony",
    max_iter       = 20,
    verbose        = FALSE
  )

  merged <- FindNeighbors(merged, reduction = "harmony", dims = HARMONY_DIMS,
                          verbose = FALSE)
  merged <- FindClusters(merged, resolution = RESOLUTION, verbose = FALSE)
  merged <- RunUMAP(merged, reduction = "harmony", dims = HARMONY_DIMS,
                    reduction.name = "umap", verbose = FALSE)
}


################################################################################
# STEP 6: Visualize
################################################################################

message("Generating UMAP plots...")

p1 <- DimPlot(merged, group.by = "sample_id",      reduction = "umap",
              label = FALSE, raster = TRUE) + ggtitle("By Sample (batch)")
p2 <- DimPlot(merged, group.by = "condition",       reduction = "umap",
              label = FALSE, raster = TRUE) + ggtitle("By Condition")
p3 <- DimPlot(merged, group.by = "technology",      reduction = "umap",
              label = FALSE, raster = TRUE) + ggtitle("By Technology")
p4 <- DimPlot(merged, group.by = "seurat_clusters", reduction = "umap",
              label = TRUE,  raster = TRUE) + ggtitle("Clusters")
p5 <- DimPlot(merged, group.by = "paper",           reduction = "umap",
              label = FALSE, raster = TRUE) + ggtitle("By Paper")
p6 <- DimPlot(merged, group.by = "gsm_id",          reduction = "umap",
              label = FALSE, raster = TRUE) + ggtitle("By GSM/GSE ID")

pdf(file.path(OUT_DIR, "UMAP_overview.pdf"), width = 18, height = 24)
print(p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 2))
dev.off()

# Known cell type annotations (from KidneyHealthy1/KidneyUUO1, EmbryosE9_5ToE13_5, and RDS atlases)
if ("cell_type_original" %in% colnames(merged@meta.data)) {
  p7 <- DimPlot(merged, group.by = "cell_type_original", reduction = "umap",
                label = TRUE, repel = TRUE, raster = TRUE) +
        ggtitle("Known cell type annotations")
  pdf(file.path(OUT_DIR, "UMAP_known_annotations.pdf"), width = 12, height = 10)
  print(p7)
  dev.off()
}


################################################################################
# STEP 7: Harmony convergence check
################################################################################

harmony_obj <- merged[["harmony"]]   # Reductions() returns names only; [[ ]] gets the object
conv_plot   <- tryCatch(harmony_obj@misc$convergence_plot,
                        error = function(e) NULL)
if (!is.null(harmony_obj) && !is.null(conv_plot)) {
  pdf(file.path(OUT_DIR, "harmony_convergence.pdf"), width = 6, height = 4)
  print(conv_plot)
  dev.off()
}


################################################################################
# STEP 8: Save
################################################################################

out_path <- file.path(OUT_DIR, "merged_harmony_integrated.rds")
message(sprintf("Saving to %s ...", out_path))
saveRDS(merged, out_path)

message("\n===== Integration complete =====")
message(sprintf("Total cells: %s", format(ncol(merged), big.mark = ",")))
message(sprintf("Clusters:    %d", length(unique(merged$seurat_clusters))))
message(sprintf("Batches:     %d samples", length(unique(merged$sample_id))))
message(sprintf("Output:      %s", OUT_DIR))
