################################################################################
# 01b_integrate_KidneyBladder_harmony.R
#
# Integrate kidney and bladder scRNA-seq urothelium cells with Harmony.
# Input files (both are scVI outputs — log-normalised expression in counts layer):
#   RenalUrothelium_uro_cells_fullgene_scvi.rds  —  11,105 kidney uro cells
#   BladderUrothelium_uro_cells_scvi.rds         —  83,715 bladder uro cells
#
# Steps:
#   1. Load & set data = counts (scVI log-normalised)
#   2. Merge (Seurat v5 per-sample layers)
#   3. JoinLayers + FindVariableFeatures (3,000 HVGs)
#   4. ScaleData (regress pct_mt)
#   5. RunPCA (30 PCs)
#   6. RunHarmony (sample_id + technology)
#   7. FindNeighbors + FindClusters (res = 0.5)
#   8. RunUMAP
#   9. Visualise → PDFs
#  10. Save KidneyBladderUrothelium_harmony_integrated.rds
################################################################################

library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)
library(dplyr)

options(future.globals.maxSize = 8 * 1024^3)

log_mem <- function(label) {
  g  <- gc(verbose = FALSE)
  mb <- sum(g[, 2])
  message(sprintf("  [mem] %s: %.1f GB in use", label, mb / 1024))
}

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/FinalUrotheliumCells"
SCR_DIR  <- file.path(BASE_DIR, "UrotheliumScripts")
OUT_DIR  <- file.path(SCR_DIR, "output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OUT_PATH <- file.path(OUT_DIR, "KidneyBladderUrothelium_harmony_integrated.rds")

# ── Parameters ────────────────────────────────────────────────────────────────
N_HVG        <- 3000
N_PCS        <- 30
HARMONY_DIMS <- 1:30
RESOLUTION   <- 0.5

# ── Input files (both scVI: counts layer = log-normalised expression) ─────────
INPUT_FILES <- list(
  list(
    path      = file.path(BASE_DIR, "RenalUrothelium_uro_cells_fullgene_scvi_metacorrection.rds"),
    sample    = "RenalUrothelium",
    tissue_lbl = "kidney"
  ),
  list(
    path      = file.path(BASE_DIR, "BladderUrothelium_uro_cells_scvi_metacorrection.rds"),
    sample    = "BladderUrothelium",
    tissue_lbl = "bladder"
  )
)



# Metadata columns expected in every object before merge
META_COLS <- c("FinalSampleId", "FinalConditionL1","FinalConditionL2","Finalscsn",  "Finaltissue", "Finalpaper", "Finalgsm_id","Finalscsn")


################################################################################
# STEP 1: Load and set data = counts (scVI log-norm)
################################################################################

message(sprintf("Loading %d scVI input files ...", length(INPUT_FILES)))

seurat_list <- lapply(INPUT_FILES, function(entry) {
  sid   <- entry$sample
  path  <- entry$path

  message(sprintf("  [scvi] %s  (%s)", sid, basename(path)))
  if (!file.exists(path)) stop("File not found: ", path)

  so <- readRDS(path)
  message(sprintf("    Loaded: %d cells × %d genes", ncol(so), nrow(so)))
  log_mem(sprintf("after loading %s", sid))

  # ── Ensure RNA assay ───────────────────────────────────────────────────────
  if (!"RNA" %in% Assays(so)) {
    cur <- DefaultAssay(so)
    message(sprintf("    Renaming assay '%s' → 'RNA'", cur))
    so  <- RenameAssays(so, assay.name = cur, new.assay.name = "RNA")
  }
  DefaultAssay(so) <- "RNA"

  # Convert to Seurat v5 Assay5 if loaded from v3/v4 format (required for JoinLayers)
  if (!inherits(so[["RNA"]], "Assay5")) {
    message("    Converting RNA assay to Seurat v5 Assay5 format")
    so[["RNA"]] <- as(so[["RNA"]], "Assay5")
  }

  # Drop non-RNA assays
  for (assay in c("SCT", "ATAC", "Spatial", "sketch", "originalexp")) {
    if (assay %in% Assays(so)) so[[assay]] <- NULL
  }

  # scVI: counts layer already holds log-normalised expression
  message("    Setting data = counts (log-norm from scVI)")
  so[["RNA"]]$data <- so[["RNA"]]$counts

  # Drop any pre-existing scale.data layers
  scale_lyrs <- grep("^scale\\.data", Layers(so), value = TRUE)
  if (length(scale_lyrs) > 0) {
    message(sprintf("    Dropping %d pre-existing scale.data layer(s)", length(scale_lyrs)))
    for (lyr in scale_lyrs) so[["RNA"]][[lyr]] <- NULL
  }

  # ── Ensure required metadata columns exist ────────────────────────────────
  for (col in META_COLS) {
    if (!col %in% colnames(so@meta.data)) {
      so@meta.data[[col]] <- NA_character_
      message(sprintf("    Added placeholder column: %s", col))
    }
  }

  # ── tissue fill-in (in case metadata is missing for some cells) ───────────
  if (any(is.na(so@meta.data$tissue))) {
    so@meta.data$tissue[is.na(so@meta.data$tissue)] <- entry$tissue_lbl
  }

  # source column (scvi for both)
  so@meta.data$source <- "scvi"

  # ── Prefix barcodes to avoid collisions after merge ───────────────────────
  so <- RenameCells(so, add.cell.id = sid)

  so
})

names(seurat_list) <- vapply(INPUT_FILES, `[[`, character(1), "sample")

message(sprintf("\nCells per input:"))
for (nm in names(seurat_list))
  message(sprintf("  %-40s : %s cells", nm,
                  format(ncol(seurat_list[[nm]]), big.mark = ",")))
message(sprintf("  %-40s : %s cells total",
                "ALL",
                format(sum(sapply(seurat_list, ncol)), big.mark = ",")))

# check the meta data for renal
seurat_list[[1]] @meta.data %>% select(all_of(META_COLS)) %>% head() %>% print()


################################################################################
# STEP 2: Merge
################################################################################

message("\nMerging objects (Seurat v5 layer-aware merge) ...")
merged <- merge(
  x          = seurat_list[[1]],
  y          = seurat_list[-1],
  merge.data = TRUE
)
rm(seurat_list); gc()

message(sprintf("Merged: %s cells × %s genes",
                format(ncol(merged), big.mark = ","),
                format(nrow(merged), big.mark = ",")))
message(sprintf("Layers: %s", paste(Layers(merged), collapse = ", ")))

# ── pct_mt ────────────────────────────────────────────────────────────────────
if (!"pct_mt" %in% colnames(merged@meta.data)) {
  message("Computing pct_mt ...")
  merged[["pct_mt"]] <- PercentageFeatureSet(merged, pattern = "^mt-|^MT-")
} else {
  n_na <- sum(is.na(merged$pct_mt))
  if (n_na > 0) {
    med_val <- median(merged$pct_mt, na.rm = TRUE)
    merged$pct_mt[is.na(merged$pct_mt)] <- med_val
    message(sprintf("Imputed %d NAs in pct_mt with median (%.4f)", n_na, med_val))
  }
}

message("\nCells per sample:")
print(sort(table(merged$sample_id), decreasing = TRUE))
message("\nCells per tissue:")
print(table(merged$tissue))
message("\nCells per condition:")
print(sort(table(merged$condition), decreasing = TRUE))


################################################################################
# STEP 3: JoinLayers + FindVariableFeatures
################################################################################

message("\nJoining layers ...")
merged <- JoinLayers(merged)

message(sprintf("FindVariableFeatures (%d HVGs) ...", N_HVG))
merged <- FindVariableFeatures(merged, selection.method = "vst",
                               nfeatures = N_HVG, verbose = FALSE)
message(sprintf("  Top 10 HVGs: %s",
                paste(head(VariableFeatures(merged), 10), collapse = ", ")))


################################################################################
# STEP 4: ScaleData
################################################################################

n_na <- sum(is.na(merged$pct_mt))
if (n_na > 0)
  merged$pct_mt[is.na(merged$pct_mt)] <- median(merged$pct_mt, na.rm = TRUE)

message("ScaleData (regress pct_mt, HVGs only) ...")
merged <- ScaleData(merged,
                    features        = VariableFeatures(merged),
                    vars.to.regress = "pct_mt",
                    verbose         = FALSE)
log_mem("after ScaleData")


################################################################################
# STEP 5: RunPCA
################################################################################

message(sprintf("RunPCA (%d PCs) ...", N_PCS))
merged <- RunPCA(merged, npcs = N_PCS, verbose = FALSE)
log_mem("after PCA")

# merged[["RNA"]]$scale.data <- NULL
gc()
log_mem("after freeing scale.data")

saveRDS(merged, file.path(OUT_DIR, "KidneyBladderUrothelium_preHarmony_PCA.rds"))
message("  Saved: KidneyBladderUrothelium_preHarmony_PCA.rds")

merged @meta.data %>% pull(sample_id) %>% table() %>% sort(decreasing = TRUE) %>% print()
################################################################################
# STEP 6: RunHarmony
################################################################################

# sample_id always; add technology if >1 level present
harmony_vars <- "sample_id"
for (v in c("Finaltechnology", "source")) {
  if (v %in% colnames(merged@meta.data) &&
      length(unique(merged@meta.data[[v]])) > 1) {
    harmony_vars <- c(harmony_vars, v)
  }
}
message(sprintf("RunHarmony (batch = %s) ...", paste(harmony_vars, collapse = " + ")))

merged <- RunHarmony(
  merged,
  group.by.vars  = harmony_vars,
  reduction      = "pca",
  reduction.save = "harmony",
  #project.dim    = FALSE,
  verbose        = FALSE
)
#merged[["pca"]] <- NULL
gc()
log_mem("after Harmony")


################################################################################
# STEP 7: FindNeighbors + FindClusters
################################################################################

message("FindNeighbors (annoy, k = 20) ...")
merged <- FindNeighbors(
  merged,
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
merged <- FindClusters(merged, resolution = RESOLUTION, verbose = FALSE)
message(sprintf("  Clusters: %d", length(unique(merged$seurat_clusters))))
gc()


################################################################################
# STEP 8: RunUMAP
################################################################################

message("RunUMAP on Harmony embedding ...")
merged <- RunUMAP(merged, reduction = "harmony", dims = HARMONY_DIMS,
                  reduction.name = "umap_harmony", verbose = FALSE)
message("  UMAP done")

saveRDS(merged, file.path(OUT_DIR, "KidneyBladderUrothelium_postHarmony_UMAP.rds"))
message("  Saved: KidneyBladderUrothelium_postHarmony_UMAP.rds")


################################################################################
# STEP 9: Visualise
################################################################################

message("Generating UMAP plots ...")

make_plot <- function(grp, title, label = FALSE, max_legend = 20) {
  if (!grp %in% colnames(merged@meta.data)) return(NULL)
  n_grp <- length(unique(merged@meta.data[[grp]]))
  p <- DimPlot(merged, group.by = grp, reduction = "umap_harmony",
               label = label, repel = label, raster = TRUE) +
    ggtitle(title)
  if (n_grp > max_legend) p <- p + NoLegend()
  p
}

plots <- Filter(Negate(is.null), list(
  make_plot("seurat_clusters", "Clusters",      label = TRUE),
  make_plot("tissue",          "By Tissue"),
  make_plot("FinalConditionL1",       "By ConditionL1"),
  make_plot("FinalConditionL2",       "By ConditionL2"),
  make_plot("FinalSampleId",       "By Sample"),
  make_plot("Finalscsn",      "By snRNA vs scRNA"),
  make_plot("Finalpaper",           "By Paper"),
  make_plot("Finalgsm_id",           "By GSM ID")
))

n_cols <- min(2, length(plots))
tryCatch({
  pdf(file.path(OUT_DIR, "KidneyBladderUrothelium_UMAP_overview.pdf"),
      width = 18, height = ceiling(length(plots) / n_cols) * 7)
  print(wrap_plots(plots, ncol = n_cols))
  dev.off()
  message("  Saved: KidneyBladderUrothelium_UMAP_overview.pdf")
}, error = function(e) {
  try(dev.off(), silent = TRUE)
  warning("Overview PDF failed: ", conditionMessage(e))
})

# Split by tissue
tryCatch({
  p_tissue <- DimPlot(merged, group.by = "seurat_clusters",
                      split.by = "tissue", reduction = "umap_harmony",
                      label = TRUE, repel = TRUE, raster = TRUE, ncol = 2) +
    ggtitle("Clusters split by tissue")
  pdf(file.path(OUT_DIR, "KidneyBladderUrothelium_UMAP_splitByTissue.pdf"),
      width = 14, height = 7)
  print(p_tissue)
  dev.off()
  message("  Saved: KidneyBladderUrothelium_UMAP_splitByTissue.pdf")
}, error = function(e) {
  try(dev.off(), silent = TRUE)
  warning("Split-by-tissue PDF failed: ", conditionMessage(e))
})

# Split by Sample
tryCatch({
  n_sample <- length(unique(merged$FinalSampleId))
  p_sample <- DimPlot(merged, group.by = "seurat_clusters",
                    split.by = "FinalSampleId", reduction = "umap_harmony",
                    label = TRUE, repel = TRUE, raster = TRUE,
                    ncol = min(4, n_sample)) +
    ggtitle("Clusters split by sample")
  pdf(file.path(OUT_DIR, "KidneyBladderUrothelium_UMAP_splitBysample.pdf"),
      width = min(4, n_sample) * 7, height = ceiling(n_sample / 4) * 7)
  print(p_sample)
  dev.off()
  message("  Saved: KidneyBladderUrothelium_UMAP_splitBysample.pdf")
}, error = function(e) {
  try(dev.off(), silent = TRUE)
  warning("Split-by-sample PDF failed: ", conditionMessage(e))
})

# Split by single cell vs single nucleus
tryCatch({
  n_tech <- length(unique(merged$Finalscsn))
  p_tech <- DimPlot(merged, group.by = "seurat_clusters",
                    split.by = "Finalscsn", reduction = "umap_harmony",
                    label = TRUE, repel = TRUE, raster = TRUE,
                    ncol = min(4, n_tech)) +
    ggtitle("Clusters split by sc vs sn")
  pdf(file.path(OUT_DIR, "KidneyBladderUrothelium_UMAP_splitByscsnRNAseq.pdf"),
      width = min(4, n_tech) * 7, height = ceiling(n_tech / 4) * 7)
  print(p_tech)
  dev.off()
  message("  Saved: KidneyBladderUrothelium_UMAP_splitByscsnRNAseq.pdf")
}, error = function(e) {
  try(dev.off(), silent = TRUE)
  warning("Split-by-scvssn PDF failed: ", technologyMessage(e))
})

# Split by technology

tryCatch({
  n_tech <- length(unique(merged$Finaltechnology))
  p_tech <- DimPlot(merged, group.by = "seurat_clusters",
                    split.by = "Finaltechnology", reduction = "umap_harmony",
                    label = TRUE, repel = TRUE, raster = TRUE,
                    ncol = min(4, n_tech)) +
    ggtitle("Clusters split by technology")
  pdf(file.path(OUT_DIR, "KidneyBladderUrothelium_UMAP_splitBytechonology.pdf"),
      width = min(4, n_tech) * 7, height = ceiling(n_tech / 4) * 7)
  print(p_tech)
  dev.off()
  message("  Saved: KidneyBladderUrothelium_UMAP_splitBytechonology.pdf")
}, error = function(e) {
  try(dev.off(), silent = TRUE)
  warning("Split-by-technology PDF failed: ", technologyMessage(e))
})


################################################################################
# STEP 10: Save
################################################################################

message(sprintf("Saving final object to %s ...", OUT_PATH))
saveRDS(merged, OUT_PATH)

message("\n===== KidneyBladder Harmony integration complete =====")
message(sprintf("Total cells : %s", format(ncol(merged), big.mark = ",")))
message(sprintf("Clusters    : %d", length(unique(merged$seurat_clusters))))
message(sprintf("Samples     : %d", length(unique(merged$FinalSampleId))))
message(sprintf("Output      : %s", OUT_PATH))


################################################################################
# STEP 11: Per-modality integration (sc vs sn)
#
# Splits the merged object by modality (single-cell vs single-nucleus),
# then runs FindVariableFeatures → ScaleData → PCA → Harmony →
# FindNeighbors/Clusters → UMAP → Visualise → Save independently for each.
#
# Modality is inferred from the 'technology' column.  Adjust SN_PATTERN below
# if your technology values differ (e.g. "snRNA-seq", "10x_nucleus", etc.).
################################################################################

SN_PATTERN <- "snRNA-seq"

merged$modality <- ifelse(
  grepl(SN_PATTERN, merged$Finalscsn, ignore.case = TRUE),
  "sn", "sc"
)
message("\nCells per modality:")
print(table(merged$modality, merged$tissue))

run_modality_pipeline <- function(so, mod_label) {
  pfx <- sprintf("%s_KidneyBladderUrothelium", mod_label)
  message(sprintf("\n========== Modality: %s  (%s cells) ==========",
                  mod_label, format(ncol(so), big.mark = ",")))
  message("  Tissues   : ", paste(sort(unique(so$tissue)),    collapse = ", "))
  message("  Samples   : ", paste(sort(unique(so$sample_id)), collapse = ", "))
  message("  Technology: ", paste(sort(unique(so$Finaltechnology)), collapse = ", "))

  # FindVariableFeatures
  message(sprintf("[%s] FindVariableFeatures (%d HVGs) ...", mod_label, N_HVG))
  so <- FindVariableFeatures(so, selection.method = "vst",
                             nfeatures = N_HVG, verbose = FALSE)

  # ScaleData
  n_na <- sum(is.na(so$pct_mt))
  if (n_na > 0)
    so$pct_mt[is.na(so$pct_mt)] <- median(so$pct_mt, na.rm = TRUE)
  message(sprintf("[%s] ScaleData (regress pct_mt) ...", mod_label))
  so <- ScaleData(so, features = VariableFeatures(so),
                  vars.to.regress = "pct_mt", verbose = FALSE)
  log_mem(sprintf("[%s] after ScaleData", mod_label))

  # RunPCA
  message(sprintf("[%s] RunPCA (%d PCs) ...", mod_label, N_PCS))
  so <- RunPCA(so, npcs = N_PCS, verbose = FALSE)
  log_mem(sprintf("[%s] after PCA", mod_label))
  saveRDS(so, file.path(OUT_DIR, sprintf("%s_preHarmony_PCA.rds", pfx)))
  message(sprintf("  Saved: %s_preHarmony_PCA.rds", pfx))

  # RunHarmony
  h_vars <- "sample_id"
  for (v in c("tissue", "technology", "source")) {
    if (v %in% colnames(so@meta.data) &&
        length(unique(so@meta.data[[v]])) > 1)
      h_vars <- c(h_vars, v)
  }
  message(sprintf("[%s] RunHarmony (batch = %s) ...",
                  mod_label, paste(h_vars, collapse = " + ")))
  so <- RunHarmony(so,
                   group.by.vars  = h_vars,
                   reduction      = "pca",
                   reduction.save = "harmony",
                   project.dim    = FALSE,
                   verbose        = FALSE)
  so[["pca"]] <- NULL
  so[["RNA"]]$scale.data <- NULL
  gc(); log_mem(sprintf("[%s] after Harmony", mod_label))

  # FindNeighbors + FindClusters
  message(sprintf("[%s] FindNeighbors (annoy, k = 20) ...", mod_label))
  so <- FindNeighbors(so, reduction = "harmony", dims = HARMONY_DIMS,
                      nn.method = "annoy", k.param = 20,
                      annoy.metric = "euclidean", n.trees = 50, verbose = FALSE)
  message(sprintf("[%s] FindClusters (res = %.2f) ...", mod_label, RESOLUTION))
  so <- FindClusters(so, resolution = RESOLUTION, verbose = FALSE)
  message(sprintf("  Clusters: %d", length(unique(so$seurat_clusters))))
  gc()

  # RunUMAP
  message(sprintf("[%s] RunUMAP ...", mod_label))
  so <- RunUMAP(so, reduction = "harmony", dims = HARMONY_DIMS,
                reduction.name = "umap_harmony", verbose = FALSE)
  saveRDS(so, file.path(OUT_DIR, sprintf("%s_postHarmony_UMAP.rds", pfx)))
  message(sprintf("  Saved: %s_postHarmony_UMAP.rds", pfx))

  # Visualise
  make_plot_mod <- function(grp, title, label = FALSE, max_legend = 20) {
    if (!grp %in% colnames(so@meta.data)) return(NULL)
    n_grp <- length(unique(so@meta.data[[grp]]))
    p <- DimPlot(so, group.by = grp, reduction = "umap_harmony",
                 label = label, repel = label, raster = TRUE) + ggtitle(title)
    if (n_grp > max_legend) p <- p + NoLegend()
    p
  }
  plots <- Filter(Negate(is.null), list(
    make_plot_mod("seurat_clusters", "Clusters",      label = TRUE),
    make_plot_mod("tissue",          "By Tissue"),
    make_plot_mod("condition",       "By Condition"),
    make_plot_mod("sample_id",       "By Sample"),
    make_plot_mod("technology",      "By Technology"),
    make_plot_mod("paper",           "By Paper")
  ))
  n_cols <- min(2, length(plots))
  tryCatch({
    pdf(file.path(OUT_DIR, sprintf("%s_UMAP_overview.pdf", pfx)),
        width = 18, height = ceiling(length(plots) / n_cols) * 7)
    print(wrap_plots(plots, ncol = n_cols)); dev.off()
    message(sprintf("  Saved: %s_UMAP_overview.pdf", pfx))
  }, error = function(e) { try(dev.off(), silent = TRUE); warning(conditionMessage(e)) })

  tryCatch({
    p_t <- DimPlot(so, group.by = "seurat_clusters", split.by = "tissue",
                   reduction = "umap_harmony", label = TRUE, repel = TRUE,
                   raster = TRUE, ncol = 2) + ggtitle("Clusters split by tissue")
    pdf(file.path(OUT_DIR, sprintf("%s_UMAP_splitByTissue.pdf", pfx)),
        width = 14, height = 7)
    print(p_t); dev.off()
    message(sprintf("  Saved: %s_UMAP_splitByTissue.pdf", pfx))
  }, error = function(e) { try(dev.off(), silent = TRUE); warning(conditionMessage(e)) })

  tryCatch({
    n_cond <- length(unique(so$condition))
    p_c <- DimPlot(so, group.by = "seurat_clusters", split.by = "condition",
                   reduction = "umap_harmony", label = TRUE, repel = TRUE,
                   raster = TRUE, ncol = min(4, n_cond)) +
      ggtitle("Clusters split by condition")
    pdf(file.path(OUT_DIR, sprintf("%s_UMAP_splitByCondition.pdf", pfx)),
        width = min(4, n_cond) * 7, height = ceiling(n_cond / 4) * 7)
    print(p_c); dev.off()
    message(sprintf("  Saved: %s_UMAP_splitByCondition.pdf", pfx))
  }, error = function(e) { try(dev.off(), silent = TRUE); warning(conditionMessage(e)) })

  # Save final
  out_path <- file.path(OUT_DIR, sprintf("%s_harmony_integrated.rds", pfx))
  message(sprintf("[%s] Saving final → %s ...", mod_label, basename(out_path)))
  saveRDS(so, out_path)

  message(sprintf("\n===== %s complete =====", pfx))
  message(sprintf("  Cells    : %s", format(ncol(so), big.mark = ",")))
  message(sprintf("  Clusters : %d", length(unique(so$seurat_clusters))))
  message(sprintf("  Output   : %s", out_path))
  invisible(so)
}

for (mod in c("sc", "sn")) {
  cells_mod <- colnames(merged)[merged$modality == mod]
  if (length(cells_mod) == 0) {
    message(sprintf("\nSkipping modality '%s': no cells found.", mod)); next
  }
  so_mod <- merged[, cells_mod]
  run_modality_pipeline(so_mod, mod)
  rm(so_mod); gc()
}

message("\n===== All modalities complete =====")

