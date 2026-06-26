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

merged <- readRDS(OUT_PATH)
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
  # plots <- Filter(Negate(is.null), list(
  #   make_plot_mod("seurat_clusters", "Clusters",      label = TRUE),
  #   make_plot_mod("tissue",          "By Tissue"),
  #   make_plot_mod("condition",       "By Condition"),
  #   make_plot_mod("sample_id",       "By Sample"),
  #   make_plot_mod("technology",      "By Technology"),
  #   make_plot_mod("paper",           "By Paper")
  # ))

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
  # mod <- c("sc")
  cells_mod <- colnames(merged)[merged$modality == mod]
  if (length(cells_mod) == 0) {
    message(sprintf("\nSkipping modality '%s': no cells found.", mod)); next
  }
  so_mod <- merged[, cells_mod]
  run_modality_pipeline(so_mod, mod)
  rm(so_mod); gc()
}

message("\n===== All modalities complete =====")

# I also want to compare the healthy bladder urothelium and healthy kidney urothelium, so I will subset the merged object to only include healthy samples from both tissues and then run Harmony integration on that subset. This will allow me to see how the healthy urothelium cells from the bladder and kidney cluster together or separately in the integrated space.

healthymerged <- subset(merged, subset = FinalConditionL1 == "HealthyBladder" | FinalConditionL1 == "HealthyKidney")

message("\nRunning Harmony on healthy subset ...")
healthymerged@meta.data$FinalConditionL1 %>% table()
healthymerged$modality <- healthymerged$FinalConditionL1
message("Cells per modality group:")
print(table(healthymerged$modality, useNA = "ifany"))

# print required parameter for run_modality_pipeline


# then run the pipeline 
run_modality_pipeline(healthymerged, "HealthyOnlyComparison")


health_out_path <- file.path(OUT_DIR, sprintf("HealthyOnlyComparison_KidneyBladderUrothelium_harmony_integrated.rds"))
healthymerged <- readRDS(health_out_path)
# 
healthymerged@meta.data %>% pull(FinalConditionL1) %>% table()

DimPlot(healthymerged, group.by = "modality", reduction = "umap_harmony", label = TRUE, repel = TRUE, raster = TRUE) + ggtitle("Healthy Bladder vs Healthy Kidney Urothelium")
# Identify the biomarkers that differed between the healthy bladder urothelium and healthy kidney urothelium. I will use the FindMarkers function in Seurat to identify differentially expressed genes between the two groups of cells. This will help me understand the molecular differences between the healthy urothelium in these two tissues.
HealthMarkersbetweenBladderKidney <- FindMarkers(healthymerged, ident.1 = "HealthyBladder", ident.2 = "HealthyKidney", group.by = "modality", logfc.threshold = 0.25, min.pct = 0.1, test.use = "MAST", verbose = TRUE)

# save the markers to a CSV file for further analysis
write.csv(HealthMarkersbetweenBladderKidney, file = file.path(OUT_DIR, "HealthyBladder_vs_HealthyKidney_markers.csv"), row.names = TRUE)

# save the RDS
saveRDS(HealthMarkersbetweenBladderKidney, file = file.path(OUT_DIR, "HealthyBladder_vs_HealthyKidney_markers.rds"))

# draw the top markers
library(tibble)
library(dplyr)

HealthMarkersbetweenBladderKidney %>% head(n=30)

top_markers_bladder <- HealthMarkersbetweenBladderKidney %>%
  tibble::rownames_to_column(var = "gene") %>%
  filter(avg_log2FC > 0, p_val_adj < 0.05, !grepl("^mt-", gene, ignore.case = TRUE)) %>%
  mutate(pct_diff = pct.1 - pct.2) %>%
  arrange(desc(pct_diff), desc(avg_log2FC)) %>%
  head(15)

top_markers_kidney <- HealthMarkersbetweenBladderKidney %>%
  tibble::rownames_to_column(var = "gene") %>%
  filter(avg_log2FC < 0, p_val_adj < 0.05, !grepl("^mt-", gene, ignore.case = TRUE)) %>%
  mutate(pct_diff = pct.2 - pct.1) %>%
  arrange(desc(pct_diff), avg_log2FC) %>%
  head(10)

top_genes <- c(top_markers_bladder$gene, top_markers_kidney$gene)

# ScaleData only ran on HVGs; scale top_genes so DoHeatmap can find them in scale.data
healthymerged <- ScaleData(healthymerged, features = top_markers_bladder$gene, verbose = FALSE)

# draw a heatmap of the top markers
Idents(healthymerged) <- "modality"

p_heatmap_bladder <- DoHeatmap(healthymerged, features = top_markers_bladder$gene, group.by = "modality",disp.max = 1.5) +
  ggtitle("Top Markers: HealthyBladder vs HealthyKidney")

ggsave(paste0(OUT_DIR, "/HealthyBladder_vs_HealthyKidney_top_bladdermarkers_heatmap.pdf"),
       plot = p_heatmap_bladder, width = 10, height = 8)

# file.path(OUT_DIR, sprintf("HealthyBladder_vs_HealthyKidney_top_markers_heatmap.pdf", pfx)

# we only selected the AP-1 and Atp genes 

top_ATP_markers_bladder <- HealthMarkersbetweenBladderKidney %>%
  tibble::rownames_to_column(var = "gene") %>%
  filter(avg_log2FC > 0, p_val_adj < 0.05, !grepl("^mt-", gene, ignore.case = TRUE)) %>%
  mutate(pct_diff = pct.1 - pct.2) %>%
  arrange(desc(pct_diff), desc(avg_log2FC)) %>% filter(
grepl("Atp", gene, ignore.case = TRUE) & pct.1 >0.8 & pct.2 < 0.2)

healthymerged <- ScaleData(healthymerged, features = top_ATP_markers_bladder$gene, verbose = FALSE)


p_ATPheatmap_bladder <- DoHeatmap(healthymerged, features = top_ATP_markers_bladder$gene, group.by = "modality",disp.max = 1.5, disp.min= -1 ) +
  ggtitle("Top Markers: HealthyBladder vs HealthyKidney")

ggsave(paste0(OUT_DIR, "/HealthyBladder_vs_HealthyKidney_top_bladderATPmarkers_heatmap.pdf"),
       plot = p_ATPheatmap_bladder, width = 10, height = 8)

# generate a dotplot for the top ATP markers, legend change to AverExp, and PerExp
p_ATPDotPlot_bladder <- DotPlot(healthymerged, features = top_ATP_markers_bladder$gene, group.by = "modality", cols = c("Spectral")) +
  ggtitle("DotPlot: HealthyBladder vs HealthyKidney top ATP markers") +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(
    color = guide_colorbar(title = "AveExp"),
    size  = guide_legend(title = "PerExp")
  )


ggsave(paste0(OUT_DIR, "/HealthyBladder_vs_HealthyKidney_top_bladderATPmarkers_Dotplot.pdf"),
       plot = p_ATPDotPlot_bladder, width = 7, height = 3.5)

# AP-1 markers
AP1Markers<- c("Fos", "Jun", "Junb", "Fosl2", "Atf3", "Egr1")

p_AP1DotPlot_bladder <- DotPlot(healthymerged, features = AP1Markers, group.by = "modality", cols = c("Spectral")) +
  ggtitle("DotPlot: HealthyBladder vs HealthyKidney top AP-1 markers") +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(
    color = guide_colorbar(title = "AveExp"),
    size  = guide_legend(title = "PerExp")
  )


ggsave(paste0(OUT_DIR, "/HealthyBladder_vs_HealthyKidney_top_bladderAP1markers_Dotplot.pdf"),
       plot = p_AP1DotPlot_bladder, width = 7, height = 3.5)


# Urothelium markers
UrothelialMarkers<- c("Krt5", "Krt14", "Krt20",  "Trp63", "Upk1a","Upk1b", "Upk2", "Upk3a","Upk3b", "Foxa1", "Gata3", "Pparg", "Krt8", "Krt18", "Krt19")

# Urothelium markers:
#   Krt8, Krt18, Krt19                         — pan-urothelium keratins
#   Upk1a, Upk1b, Upk2, Upk3a, Upk3b          — uroplakins
#   Krt20, Krt5, Krt14, Trp63                  — umbrella / basal markers
#   Foxa1, Gata3, Pparg                        — urothelial TFs

p_UrothelialDotPlot_bladder <- DotPlot(healthymerged, features = UrothelialMarkers, group.by = "modality", cols = c("Spectral")) +
  ggtitle("DotPlot: HealthyBladder vs HealthyKidney Urothelium markers") +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(
    color = guide_colorbar(title = "AveExp"),
    size  = guide_legend(title = "PerExp")
  )


ggsave(paste0(OUT_DIR, "/HealthyBladder_vs_HealthyKidney_Urotheliummarkers_Dotplot.pdf"),
       plot = p_UrothelialDotPlot_bladder, width = 8, height = 3.5)


CandiategenesMarkers <- c(
# Renal identity
"Pax8",
"Pax2",
"Glis3",
"Fgfr2",

# Cilia / collecting system
"Pkhd1",
"Bicc1",

# Epithelial architecture
"Magi1",
"Cgnl1",
"Ptpn14",

# Basement membrane
"Col4a3",
"Col4a4",
"Col4a5",

# Bladder AP-1 program
"Fos",
"Jun",
"Junb",
"Fosl2", 
"Atf3",

# Bladder metabolic program
"Atp5g2",
"Atp5l",
"Atpif1",
"Atp5e"
)

p_CandidateDotPlot_bladderkidney <- DotPlot(healthymerged, features = CandiategenesMarkers, group.by = "modality", cols = c("Spectral")) +
  ggtitle("DotPlot: HealthyBladder vs HealthyKidney Candidate markers") +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(
    color = guide_colorbar(title = "AveExp"),
    size  = guide_legend(title = "PerExp")
  )


ggsave(paste0(OUT_DIR, "/HealthyBladder_vs_HealthyKidney_FinalCandidateMmarkers_Dotplot.pdf"),
       plot = p_CandidateDotPlot_bladderkidney, width = 10, height = 3.5)

CandiategenesKidneyMarkers <- c(
# Renal identity
"Pax8",
"Pax2",
"Glis3",
"Fgfr2",

# Cilia / collecting system
"Pkhd1",
"Bicc1",

# Epithelial architecture
"Magi1",
"Cgnl1",
"Ptpn14",

# Basement membrane
"Col4a3",
"Col4a4",
"Col4a5") 
p_CandidateDotPlot_kidney <- DotPlot(healthymerged, features = CandiategenesKidneyMarkers, group.by = "modality", cols = c("Spectral")) +
  # ggtitle("DotPlot: HealthyBladder vs HealthyKidney Candidate markers") +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(
    color = guide_colorbar(title = "AveExp"),
    size  = guide_legend(title = "PerExp")
  )


ggsave(paste0(OUT_DIR, "/HealthyBladder_vs_HealthyKidney_FinalCandidateKidneyMmarkers_Dotplot.pdf"),
       plot = p_CandidateDotPlot_kidney, width = 8, height = 4)

# I also want to put the kidney epithelial and immune cells markers for healthy kidney and healthy bladder together to see how they cluster. I will subset the merged object to only include epithelial and immune cells from both tissues, then run Harmony integration on that subset. This will allow me to see how the epithelial and immune cells from the healthy bladder and healthy kidney cluster together or separately in the integrated space.

# tubule cells markers
RenalEpimuneMarkers <- c(
  "Slc34a1", "Lrp2",    "Cubn",     # proximal tubule
  "Umod",    "Slc12a1",             # TAL
  "Slc12a3", "Pvalb",               # DCT
  "Calb1", "Trpv5", "Atp2b4",       # CNT
  "Aqp2",     "Scnn1g",   # collecting duct principal
  "Atp6v1b1", "Slc4a1",  "Foxi1",   # intercalated cells
  "Ptprc", "Cd3e", "Cd4", "Cd8a", "Cd14", "Cd68", "Itgam"   # Immune markers
)

p_tubuleImmDotPlot_kidney <- DotPlot(healthymerged, features = RenalEpimuneMarkers, group.by = "modality", cols = c("Spectral")) +
  # ggtitle("DotPlot: HealthyBladder vs HealthyKidney Candidate markers") +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(
    color = guide_colorbar(title = "AveExp"),
    size  = guide_legend(title = "PerExp")
  )


ggsave(paste0(OUT_DIR, "/HealthyBladder_vs_HealthyKidney_FinalTubuleImmuneKidneyMmarkers_Dotplot.pdf"),
       plot = p_tubuleImmDotPlot_kidney, width = 8, height = 4)


################################################################################
# STEP 12: GO & KEGG enrichment — HealthyBladder vs HealthyKidney urothelium
#
# DEG filters (from HealthMarkersbetweenBladderKidney):
#   - exclude mt genes
#   - FDR (p_val_adj) < 0.05
#   - |avg_log2FC| > 0.25
#   - expressed in ≥10% of cells in at least one group (pct.1 ≥ 0.1 | pct.2 ≥ 0.1)
#
# Gene sets:
#   bladder_up : avg_log2FC > 0.25  (higher in bladder)
#   kidney_up  : avg_log2FC < -0.25 (higher in kidney)
################################################################################

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)

message("\n===== STEP 12: GO & KEGG enrichment =====")

# ── Filter DEGs ───────────────────────────────────────────────────────────────
# Exclude mt genes (^mt-) and ribosomal protein genes (^Rps, ^Rpl, ^Mrps, ^Mrpl)
is_excluded_gene <- function(g) {
  grepl("^mt-", g, ignore.case = TRUE) |
  grepl("^Rps", g, ignore.case = TRUE) |
  grepl("^Rpl", g, ignore.case = TRUE) |
  grepl("^Mrps", g, ignore.case = TRUE) |
  grepl("^Mrpl", g, ignore.case = TRUE)
}

deg_all <- HealthMarkersbetweenBladderKidney %>%
  tibble::rownames_to_column("gene") %>%
  filter(
    !is_excluded_gene(gene),
    p_val_adj < 0.05,
    abs(avg_log2FC) > 0.25,
    pct.1 >= 0.1 | pct.2 >= 0.1
  )

message(sprintf("  DEGs after filtering: %d total (%d bladder-up, %d kidney-up)",
  nrow(deg_all),
  sum(deg_all$avg_log2FC > 0),
  sum(deg_all$avg_log2FC < 0)))

bladder_genes <- deg_all %>% filter(avg_log2FC >  0.25) %>% pull(gene)
kidney_genes  <- deg_all %>% filter(avg_log2FC < -0.25) %>% pull(gene)

# Background = all genes tested in FindMarkers (excluding mt and ribosomal)
background_genes <- HealthMarkersbetweenBladderKidney %>%
  tibble::rownames_to_column("gene") %>%
  filter(!is_excluded_gene(gene)) %>%
  pull(gene)

# ── Symbol → Entrez ID conversion ────────────────────────────────────────────
sym2entrez <- function(genes) {
  mapped <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Mm.eg.db, drop = TRUE)
  mapped$ENTREZID
}

bladder_entrez    <- sym2entrez(bladder_genes)
kidney_entrez     <- sym2entrez(kidney_genes)
background_entrez <- sym2entrez(background_genes)

message(sprintf("  Bladder-up: %d genes → %d Entrez IDs", length(bladder_genes), length(bladder_entrez)))
message(sprintf("  Kidney-up : %d genes → %d Entrez IDs", length(kidney_genes),  length(kidney_entrez)))

# ── GO enrichment (Biological Process) ───────────────────────────────────────
run_go <- function(entrez_ids, label, universe = background_entrez) {
  if (length(entrez_ids) < 5) {
    message(sprintf("  [GO] %s: too few IDs (%d), skipping", label, length(entrez_ids)))
    return(NULL)
  }
  ego <- enrichGO(
    gene          = entrez_ids,
    universe      = universe,
    OrgDb         = org.Mm.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1,
    readable      = TRUE
  )
  message(sprintf("  [GO BP] %s: %d enriched terms", label, nrow(as.data.frame(ego))))
  ego
}

go_bladder <- run_go(bladder_entrez, "bladder_up")
go_kidney  <- run_go(kidney_entrez,  "kidney_up")

# Save GO results
if (!is.null(go_bladder) && nrow(as.data.frame(go_bladder)) > 0) {
  write.csv(as.data.frame(go_bladder),
            file.path(OUT_DIR, "GO_BP_BladderUp_vs_KidneyUrothelium.csv"),
            row.names = FALSE)

  p_go_bladder_dot <- dotplot(go_bladder, showCategory = 20, font.size = 9) +
    ggtitle("GO BP: Bladder-up genes (vs Kidney urothelium)")
  ggsave(file.path(OUT_DIR, "GO_BP_BladderUp_dotplot.pdf"),
         plot = p_go_bladder_dot, width = 9, height = 10)

  p_go_bladder_bar <- barplot(go_bladder, showCategory = 20, font.size = 9) +
    ggtitle("GO BP: Bladder-up genes (vs Kidney urothelium)")
  ggsave(file.path(OUT_DIR, "GO_BP_BladderUp_barplot.pdf"),
         plot = p_go_bladder_bar, width = 9, height = 10)

  message("  Saved GO BP results for bladder-up genes")
}

if (!is.null(go_kidney) && nrow(as.data.frame(go_kidney)) > 0) {
  write.csv(as.data.frame(go_kidney),
            file.path(OUT_DIR, "GO_BP_KidneyUp_vs_BladderUrothelium.csv"),
            row.names = FALSE)

  p_go_kidney_dot <- dotplot(go_kidney, showCategory = 20, font.size = 9) +
    ggtitle("GO BP: Kidney-up genes (vs Bladder urothelium)")
  ggsave(file.path(OUT_DIR, "GO_BP_KidneyUp_dotplot.pdf"),
         plot = p_go_kidney_dot, width = 9, height = 10)

  p_go_kidney_bar <- barplot(go_kidney, showCategory = 20, font.size = 9) +
    ggtitle("GO BP: Kidney-up genes (vs Bladder urothelium)")
  ggsave(file.path(OUT_DIR, "GO_BP_KidneyUp_barplot.pdf"),
         plot = p_go_kidney_bar, width = 9, height = 10)

  message("  Saved GO BP results for kidney-up genes")
}

# Combined GO dotplot (top 15 per direction)
if (!is.null(go_bladder) && !is.null(go_kidney) &&
    nrow(as.data.frame(go_bladder)) > 0 && nrow(as.data.frame(go_kidney)) > 0) {
  go_bladder_df <- as.data.frame(go_bladder) %>%
    arrange(p.adjust) %>% head(15) %>%
    mutate(direction = "Bladder-up")
  go_kidney_df <- as.data.frame(go_kidney) %>%
    arrange(p.adjust) %>% head(15) %>%
    mutate(direction = "Kidney-up")

  go_combined <- bind_rows(go_bladder_df, go_kidney_df) %>%
    mutate(
      Description = factor(Description, levels = rev(unique(Description))),
      log10_padj  = -log10(p.adjust)
    )

  p_go_combined <- ggplot(go_combined,
    aes(x = log10_padj, y = Description, fill = direction, size = Count)) +
    geom_point(shape = 21, alpha = 0.85) +
    scale_fill_manual(values = c("Bladder-up" = "#E45F5F", "Kidney-up" = "#4F8FD6")) +
    scale_size_continuous(range = c(3, 10)) +
    labs(x = "-log10(FDR)", y = NULL,
         title = "GO BP: Bladder vs Kidney urothelium (top 15 per direction)",
         fill = "Direction", size = "Gene count") +
    theme_bw(base_size = 10) +
    theme(axis.text.y = element_text(size = 8))

  ggsave(file.path(OUT_DIR, "GO_BP_Combined_BladderKidney_dotplot.pdf"),
         plot = p_go_combined, width = 11, height = 12)
  message("  Saved combined GO BP dotplot")
}

# ── KEGG enrichment ───────────────────────────────────────────────────────────
run_kegg <- function(entrez_ids, label, universe = background_entrez) {
  if (length(entrez_ids) < 5) {
    message(sprintf("  [KEGG] %s: too few IDs (%d), skipping", label, length(entrez_ids)))
    return(NULL)
  }
  ekegg <- enrichKEGG(
    gene          = entrez_ids,
    universe      = universe,
    organism      = "mmu",
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1
  )
  message(sprintf("  [KEGG] %s: %d enriched pathways", label, nrow(as.data.frame(ekegg))))
  ekegg
}

kegg_bladder <- run_kegg(bladder_entrez, "bladder_up")
kegg_kidney  <- run_kegg(kidney_entrez,  "kidney_up")

# Save KEGG results
if (!is.null(kegg_bladder) && nrow(as.data.frame(kegg_bladder)) > 0) {
  write.csv(as.data.frame(kegg_bladder),
            file.path(OUT_DIR, "KEGG_BladderUp_vs_KidneyUrothelium.csv"),
            row.names = FALSE)

  p_kegg_bladder_dot <- dotplot(kegg_bladder, showCategory = 20, font.size = 9) +
    ggtitle("KEGG: Bladder-up genes (vs Kidney urothelium)")
  ggsave(file.path(OUT_DIR, "KEGG_BladderUp_dotplot.pdf"),
         plot = p_kegg_bladder_dot, width = 9, height = 8)

  p_kegg_bladder_bar <- barplot(kegg_bladder, showCategory = 20, font.size = 9) +
    ggtitle("KEGG: Bladder-up genes (vs Kidney urothelium)")
  ggsave(file.path(OUT_DIR, "KEGG_BladderUp_barplot.pdf"),
         plot = p_kegg_bladder_bar, width = 9, height = 8)

  message("  Saved KEGG results for bladder-up genes")
}

if (!is.null(kegg_kidney) && nrow(as.data.frame(kegg_kidney)) > 0) {
  write.csv(as.data.frame(kegg_kidney),
            file.path(OUT_DIR, "KEGG_KidneyUp_vs_BladderUrothelium.csv"),
            row.names = FALSE)

  p_kegg_kidney_dot <- dotplot(kegg_kidney, showCategory = 20, font.size = 9) +
    ggtitle("KEGG: Kidney-up genes (vs Bladder urothelium)")
  ggsave(file.path(OUT_DIR, "KEGG_KidneyUp_dotplot.pdf"),
         plot = p_kegg_kidney_dot, width = 9, height = 8)

  p_kegg_kidney_bar <- barplot(kegg_kidney, showCategory = 20, font.size = 9) +
    ggtitle("KEGG: Kidney-up genes (vs Bladder urothelium)")
  ggsave(file.path(OUT_DIR, "KEGG_KidneyUp_barplot.pdf"),
         plot = p_kegg_kidney_bar, width = 9, height = 8)

  message("  Saved KEGG results for kidney-up genes")
}

# Combined KEGG dotplot (top 15 per direction)
if (!is.null(kegg_bladder) && !is.null(kegg_kidney) &&
    nrow(as.data.frame(kegg_bladder)) > 0 && nrow(as.data.frame(kegg_kidney)) > 0) {
  kegg_bladder_df <- as.data.frame(kegg_bladder) %>%
    arrange(p.adjust) %>% head(15) %>%
    mutate(direction = "Bladder-up")
  kegg_kidney_df <- as.data.frame(kegg_kidney) %>%
    arrange(p.adjust) %>% head(15) %>%
    mutate(direction = "Kidney-up")

  kegg_combined <- bind_rows(kegg_bladder_df, kegg_kidney_df) %>%
    mutate(
      Description = factor(Description, levels = rev(unique(Description))),
      log10_padj  = -log10(p.adjust)
    )

  p_kegg_combined <- ggplot(kegg_combined,
    aes(x = log10_padj, y = Description, fill = direction, size = Count)) +
    geom_point(shape = 21, alpha = 0.85) +
    scale_fill_manual(values = c("Bladder-up" = "#E45F5F", "Kidney-up" = "#4F8FD6")) +
    scale_size_continuous(range = c(3, 10)) +
    labs(x = "-log10(FDR)", y = NULL,
         title = "KEGG: Bladder vs Kidney urothelium (top 15 per direction)",
         fill = "Direction", size = "Gene count") +
    theme_bw(base_size = 10) +
    theme(axis.text.y = element_text(size = 8))

  ggsave(file.path(OUT_DIR, "KEGG_Combined_BladderKidney_dotplot.pdf"),
         plot = p_kegg_combined, width = 11, height = 10)
  message("  Saved combined KEGG dotplot")
}

message("\n===== GO & KEGG enrichment complete =====")
message(sprintf("  Outputs in: %s", OUT_DIR))
