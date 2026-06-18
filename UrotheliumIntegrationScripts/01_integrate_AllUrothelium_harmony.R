################################################################################
# 01_integrate_AllUrothelium_harmony.R
#
# Integrate all urothelium samples with Harmony:
#   BladderUrothelium_uro_cells_scvi.rds     — bladder urothelium (scVI)
#   RenalUrothelium_uro_cells_fullgene_scvi.rds — kidney urothelium (scVI)
#   BladderHomogenate1_qc.rds               — bladder urothelium organoid
#   BladderHomogenate2_qc.rds               — bladder urothelium organoid
#   KudoUUOUrothelium_qc.rds               — kidney urothelium organoid (KUDO)
#   MouseUreterRecon1_qc.rds               — ureter organoid
#
# scVI files: counts layer = log-normalised expression → set data = counts
# _qc files:  counts layer = raw counts               → NormalizeData
#
# Steps:
#   1. Load & per-sample normalise
#   2. Harmonise metadata
#   3. Merge (Seurat v5 per-sample layers)
#   4. JoinLayers + FindVariableFeatures (3 000 HVGs)
#   5. ScaleData (regress pct_mt)
#   6. RunPCA (20 PCs)
#   7. RunHarmony (sample_id + technology)
#   8. FindNeighbors + FindClusters (res = 0.5)
#   9. RunUMAP
#  10. Visualise → PDFs
#  11. Save AllUrothelium_harmony_integrated.rds
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

OUT_PATH <- file.path(OUT_DIR, "AllUrothelium_harmony_integrated.rds")

# ── Parameters ────────────────────────────────────────────────────────────────
N_HVG        <- 3000
N_PCS        <- 20
HARMONY_DIMS <- 1:20
RESOLUTION   <- 0.5

# ── Input files ───────────────────────────────────────────────────────────────
# Two types:
#   scvi  — log-normalised X stored in counts layer; already batch-corrected
#   qc    — raw counts; requires NormalizeData
INPUT_FILES <- list(
  list(path = file.path(BASE_DIR, "BladderUrothelium_uro_cells_scvi.rds"),
       type = "scvi",  sample = "BladderUrothelium"),
  list(path = file.path(BASE_DIR, "RenalUrothelium_uro_cells_fullgene_scvi.rds"),
       type = "scvi",  sample = "RenalUrothelium"),
  list(path = file.path(BASE_DIR, "BladderHomogenate1_qc.rds"),
       type = "qc",    sample = "BladderHomogenate1"),
  list(path = file.path(BASE_DIR, "BladderHomogenate2_qc.rds"),
       type = "qc",    sample = "BladderHomogenate2"),
  list(path = file.path(BASE_DIR, "KudoUUOUrothelium_qc.rds"),
       type = "qc",    sample = "KudoUUOUrothelium"),
  list(path = file.path(BASE_DIR, "MouseUreterRecon1_qc.rds"),
       type = "qc",    sample = "MouseUreterRecon1")
)

# ── Default metadata for _qc organoid samples (add if columns are missing) ───
# Modify these values to match the correct paper/GSM IDs for each dataset.
QC_META_DEFAULTS <- data.frame(
  sample_id  = c("BladderHomogenate1", "BladderHomogenate2",
                 "KudoUUOUrothelium",  "MouseUreterRecon1"),
  condition  = c("BladderOrganoid",    "BladderOrganoid",
                 "KidneyOrganoid",     "UreterOrganoid"),
  technology = c("10X",                "10X",
                 "PIPseq",             "10X"),
  tissue     = c("bladder",            "bladder",
                 "kidney",             "ureter"),
  paper      = c("PMID_31562298",      "PMID_31562298",
                 "Kudo2026",           "Unknown"),
  gsm_id     = c("GSM3827175",         "GSM3827176",
                 "CustomedKudoUrothelium", "CustomedUreter"),
  stringsAsFactors = FALSE
)
rownames(QC_META_DEFAULTS) <- QC_META_DEFAULTS$sample_id

# Metadata columns that must be present in every object before merge
META_COLS <- c("sample_id", "condition", "technology", "tissue", "paper", "gsm_id", "source")

# ── Helper: fill in missing metadata from defaults table ─────────────────────
add_missing_meta <- function(so, sid, defaults_row) {
  for (col in setdiff(names(defaults_row), "sample_id")) {
    if (!col %in% colnames(so@meta.data) ||
        all(is.na(so@meta.data[[col]]))) {
      so@meta.data[[col]] <- defaults_row[[col]]
      message(sprintf("    Added missing metadata column '%s' = '%s'",
                      col, defaults_row[[col]]))
    }
  }
  if (!"sample_id" %in% colnames(so@meta.data) ||
      all(is.na(so@meta.data$sample_id))) {
    so@meta.data$sample_id <- sid
  }
  so
}

# ── Helper: handle KudoUUOUrothelium per-sub-sample metadata ─────────────────
# The KUDO object may have a Sample column with GOF/LOF/Vehicle/Yoda
# (organoid conditions) and/or Health/UUO_Day2/… (sci-RNA-seq conditions).
harmonise_kudo_meta <- function(so) {
  if (!"Sample" %in% colnames(so@meta.data)) return(so)
  message("  KudoUUOUrothelium: splitting sample_id by Sample column ...")

  so@meta.data <- so@meta.data %>%
    mutate(
      sample_id = case_when(
        sample_id == "KudoUUOUrothelium" ~ paste0("KUDO_", Sample),
        TRUE ~ sample_id
      ),
      technology = case_when(
        Sample %in% c("GOF", "LOF", "Vehicle", "Yoda") ~ "PIPseq",
        Sample %in% c("Health", "UUO_Day2", "UUO_Day4",
                      "UUO_Day6", "UUO_Day10", "UUO_Day14") ~ "sci-RNA-seq3",
        TRUE ~ technology
      ),
      condition = case_when(
        Sample %in% c("GOF", "LOF", "Vehicle", "Yoda") ~ "KidneyUrogenoid_Organoid",
        Sample == "Health"     ~ "Healthy",
        Sample == "UUO_Day2"  ~ "UUO_2days",
        Sample == "UUO_Day4"  ~ "UUO_4days",
        Sample == "UUO_Day6"  ~ "UUO_6days",
        Sample == "UUO_Day10" ~ "UUO_10days",
        Sample == "UUO_Day14" ~ "UUO_14days",
        TRUE ~ condition
      ),
      paper = case_when(
        Sample %in% c("GOF", "LOF", "Vehicle", "Yoda") ~ "Kudo2026",
        Sample %in% c("Health", "UUO_Day2", "UUO_Day4",
                      "UUO_Day6", "UUO_Day10", "UUO_Day14") ~ "PMID_36265491",
        TRUE ~ paper
      )
    )
  so
}


################################################################################
# STEP 1: Load and per-sample normalise
################################################################################

message(sprintf("Loading %d input files ...", length(INPUT_FILES)))

seurat_list <- lapply(INPUT_FILES, function(entry) {
  sid  <- entry$sample
  ftyp <- entry$type
  path <- entry$path

  message(sprintf("  [%s] %s  (%s)", ftyp, sid, basename(path)))
  if (!file.exists(path))
    stop("File not found: ", path)

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

  # Drop non-RNA assays to save memory
  for (assay in c("SCT", "ATAC", "Spatial", "sketch", "originalexp")) {
    if (assay %in% Assays(so)) so[[assay]] <- NULL
  }

  # ── Normalisation ─────────────────────────────────────────────────────────
  if (ftyp == "scvi") {
    # counts layer already holds log-normalised expression from scVI h5ad
    message("    scVI file: setting data = counts (log-norm from scVI)")
    so[["RNA"]]$data <- so[["RNA"]]$counts
  } else {
    # _qc file: raw counts → log-normalise
    message("    QC file: running NormalizeData ...")
    so <- NormalizeData(so, normalization.method = "LogNormalize",
                        scale.factor = 10000, verbose = FALSE)
  }

  # Drop any pre-existing scale.data layers so they don't conflict with the
  # post-merge ScaleData call (KudoUUOUrothelium ships with per-sample layers)
  scale_lyrs <- grep("^scale\\.data", Layers(so), value = TRUE)
  if (length(scale_lyrs) > 0) {
    message(sprintf("    Dropping %d pre-existing scale.data layer(s)", length(scale_lyrs)))
    for (lyr in scale_lyrs) so[["RNA"]][[lyr]] <- NULL
  }

  # ── source column ─────────────────────────────────────────────────────────
  so@meta.data$source <- ftyp

  # ── Metadata defaults for _qc organoid samples ────────────────────────────
  if (ftyp == "qc" && sid %in% rownames(QC_META_DEFAULTS)) {
    so <- add_missing_meta(so, sid, QC_META_DEFAULTS[sid, ])
  }

  # ── KudoUUOUrothelium: per-sub-sample metadata ────────────────────────────
  if (sid == "KudoUUOUrothelium") {
    so <- harmonise_kudo_meta(so)
  }

  # ── Prefix barcodes with sample name to avoid collisions ──────────────────
  so <- RenameCells(so, add.cell.id = sid)

  so
})

names(seurat_list) <- vapply(INPUT_FILES, `[[`, character(1), "sample")
seurat_list

message(sprintf("\nTotal cells loaded: %s",
                format(sum(sapply(seurat_list, ncol)), big.mark = ",")))
for (nm in names(seurat_list))
  message(sprintf("  %-40s : %s cells", nm,
                  format(ncol(seurat_list[[nm]]), big.mark = ",")))


################################################################################
# STEP 2: Merge
################################################################################

message("\nMerging all objects (Seurat v5 layer-aware merge) ...")
merged <- merge(
  x            = seurat_list[[1]],
  y            = seurat_list[-1],
  add.cell.ids = NULL,
  merge.data   = TRUE
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
if ("condition" %in% colnames(merged@meta.data)) {
  message("\nCells per condition:")
  print(sort(table(merged$condition), decreasing = TRUE))
}


################################################################################
# STEP 3: JoinLayers + FindVariableFeatures
################################################################################

message("\nJoining layers for ScaleData/PCA ...")
merged <- JoinLayers(merged)

message(sprintf("FindVariableFeatures (%d HVGs) ...", N_HVG))
merged <- FindVariableFeatures(merged, selection.method = "vst",
                               nfeatures = N_HVG, verbose = FALSE)
message(sprintf("  Top 10 HVGs: %s",
                paste(head(VariableFeatures(merged), 10), collapse = ", ")))


################################################################################
# STEP 4: ScaleData
################################################################################

# Fill any remaining pct_mt NAs before regression
n_na <- sum(is.na(merged$pct_mt))
if (n_na > 0) {
  merged$pct_mt[is.na(merged$pct_mt)] <- median(merged$pct_mt, na.rm = TRUE)
}

message("ScaleData (regress pct_mt, HVGs only) ...")
merged <- ScaleData(merged,
                    features        = VariableFeatures(merged),
                    vars.to.regress = "pct_mt",
                    verbose         = FALSE)
message("  ScaleData done")
log_mem("after ScaleData")


################################################################################
# STEP 5: RunPCA
################################################################################

message(sprintf("RunPCA (%d PCs) ...", N_PCS))
merged <- RunPCA(merged, npcs = N_PCS, verbose = FALSE)
message("  PCA done")
log_mem("after PCA")

merged[["RNA"]]$scale.data <- NULL
gc()
log_mem("after freeing scale.data")

# we will save the object with the PCA embedding
saveRDS(merged, file.path(OUT_DIR, "AllUrothelium_preHarmony_PCA.rds"))

################################################################################
# STEP 6: RunHarmony
################################################################################

# Use sample_id always; add technology if >1 level
harmony_vars <- "sample_id"
for (v in c("technology", "source")) {
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
  verbose        = FALSE
)
merged[["pca"]] <- NULL
gc()
message("  Harmony done")
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

# save the object with the UMAP embedding
saveRDS(merged, file.path(OUT_DIR, "AllUrothelium_postHarmony_UMAP.rds"))

################################################################################
# STEP 9: Visualise
################################################################################

message("Generating UMAP plots ...")

# Suppress legend for high-cardinality groups to prevent viewport overflow
make_plot <- function(grp, title, label = FALSE, max_legend = 20) {
  if (!grp %in% colnames(merged@meta.data)) return(NULL)
  n_grp <- length(unique(merged@meta.data[[grp]]))
  p <- DimPlot(merged, group.by = grp, reduction = "umap_harmony",
               label = label, repel = label, raster = TRUE) + ggtitle(title)
  if (n_grp > max_legend) p <- p + NoLegend()
  p
}

plots <- Filter(Negate(is.null), list(
  make_plot("seurat_clusters", "Clusters",      label = TRUE),
  make_plot("condition",       "By Condition"),
  make_plot("sample_id",       "By Sample"),
  make_plot("technology",      "By Technology"),
  make_plot("tissue",          "By Tissue"),
  make_plot("paper",           "By Paper"),
  make_plot("source",          "By Source (scvi vs qc)")
))

n_cols <- min(2, length(plots))
tryCatch({
  pdf(file.path(OUT_DIR, "AllUrothelium_UMAP_overview.pdf"),
      width = 18, height = ceiling(length(plots) / n_cols) * 7)
  print(wrap_plots(plots, ncol = n_cols))
  dev.off()
  message("  Saved: AllUrothelium_UMAP_overview.pdf")
}, error = function(e) {
  try(dev.off(), silent = TRUE)
  warning("Overview PDF failed: ", conditionMessage(e))
})

# UMAP split by tissue
if ("tissue" %in% colnames(merged@meta.data)) {
  tryCatch({
    p_tissue <- DimPlot(merged, group.by = "seurat_clusters",
                        split.by = "tissue", reduction = "umap_harmony",
                        label = TRUE, repel = TRUE, raster = TRUE, ncol = 3) +
      ggtitle("Clusters split by tissue")
    pdf(file.path(OUT_DIR, "AllUrothelium_UMAP_splitByTissue.pdf"),
        width = 21, height = 7)
    print(p_tissue)
    dev.off()
    message("  Saved: AllUrothelium_UMAP_splitByTissue.pdf")
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    warning("Split-by-tissue PDF failed: ", conditionMessage(e))
  })
}


################################################################################
# STEP 10: Save
################################################################################

message(sprintf("Saving to %s ...", OUT_PATH))
saveRDS(merged, OUT_PATH)

message("\n===== AllUrothelium Harmony integration complete =====")
message(sprintf("Total cells : %s", format(ncol(merged), big.mark = ",")))
message(sprintf("Clusters    : %d", length(unique(merged$seurat_clusters))))
message(sprintf("Samples     : %d", length(unique(merged$sample_id))))
message(sprintf("Output      : %s", OUT_PATH))
