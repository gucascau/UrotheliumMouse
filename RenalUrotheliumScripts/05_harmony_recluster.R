################################################################################
# 05_harmony_recluster.R — Cluster renal urothelial cells using scVI latent space
#
# Input : output/RenalUrothelium_uro_cells_fullgene_scvi.rds
#         counts layer = raw counts
#         data layer   = log-normalised expression from AnnData layers["lognorm"]
#         X_scVI reduction = scVI latent space (already batch-corrected)
#
# Steps:
#   1. Load Seurat object
#   2. Verify counts/data layers + X_scVI reduction
#   3. Compute pct_mt from raw counts
#   4. FindNeighbors + FindClusters (res = 0.5) on X_scVI
#   5. RunUMAP on X_scVI
#   6. UMAP + marker plots → PDFs
#   7. Save output/RenalUrothelium_uro_cells_fullgene_harmony_integrated.rds
################################################################################

library(Seurat)
library(ggplot2)
library(patchwork)

options(future.globals.maxSize = 4 * 1024^3)

log_mem <- function(label) {
  g  <- gc(verbose = FALSE)
  mb <- sum(g[, 2])
  message(sprintf("  [mem] %s: %.1f GB in use", label, mb / 1024))
}

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
OUT_DIR  <- file.path(BASE_DIR, "RenalUrotheliumScripts", "output")
dir.create(OUT_DIR, showWarnings = FALSE)

IN_PATH  <- file.path(OUT_DIR, "RenalUrothelium_uro_cells_fullgene_scvi.rds")
OUT_PATH <- file.path(OUT_DIR, "RenalUrothelium_uro_cells_fullgene_harmony_integrated.rds")

# ── Parameters ────────────────────────────────────────────────────────────────
RESOLUTION <- 0.5


################################################################################
# STEP 1: Load
################################################################################

if (!file.exists(IN_PATH))
  stop("Input not found: ", IN_PATH)

message("Loading RenalUrothelium_uro_cells_fullgene_scvi.rds ...")
so <- readRDS(IN_PATH)
message(sprintf("  Loaded: %d cells × %d genes", ncol(so), nrow(so)))
log_mem("after load")

for (col in c("sample_id", "condition", "technology")) {
  if (col %in% colnames(so@meta.data)) {
    message(sprintf("  Cells per %s:", col))
    print(sort(table(so@meta.data[[col]]), decreasing = TRUE))
  }
}


################################################################################
# STEP 2: Verify layers + X_scVI reduction
################################################################################

message("Verifying RNA assay layers ...")
rna_layers <- Layers(so[["RNA"]])
message("  RNA layers: ", paste(rna_layers, collapse = ", "))

for (layer in c("counts", "data")) {
  if (!layer %in% rna_layers)
    stop("RNA assay is missing required layer: ", layer)
}

message("Available reductions: ", paste(Reductions(so), collapse = ", "))
if (!"X_scVI" %in% Reductions(so))
  stop("X_scVI reduction not found. Available: ", paste(Reductions(so), collapse = ", "))

SCVI_DIMS <- seq_len(ncol(Embeddings(so, "X_scVI")))
message(sprintf("  X_scVI dims: %d", length(SCVI_DIMS)))
log_mem("after layer check")


################################################################################
# STEP 3: pct_mt
################################################################################

message("Computing pct_mt ...")
so[["pct_mt"]] <- PercentageFeatureSet(so, pattern = "^mt-")
message(sprintf("  pct_mt range: %.2f – %.2f", min(so$pct_mt), max(so$pct_mt)))

n_na <- sum(is.na(so$pct_mt))
if (n_na > 0) {
  med_val <- median(so$pct_mt, na.rm = TRUE)
  so$pct_mt[is.na(so$pct_mt)] <- med_val
  message(sprintf("  Imputed %d NAs in pct_mt with median (%.4f)", n_na, med_val))
}


################################################################################
# STEP 4: FindNeighbors + FindClusters on X_scVI
################################################################################

message("FindNeighbors on X_scVI (annoy, k=20) ...")
so <- FindNeighbors(
  so,
  reduction    = "X_scVI",
  dims         = SCVI_DIMS,
  nn.method    = "annoy",
  k.param      = 20,
  annoy.metric = "euclidean",
  n.trees      = 50,
  verbose      = FALSE
)
log_mem("after FindNeighbors")

message(sprintf("FindClusters (resolution = %.2f) ...", RESOLUTION))
snn_graphs <- grep("_snn$", names(so@graphs), value = TRUE)
if (length(snn_graphs) == 0)
  stop("No SNN graph found. Available graphs: ", paste(names(so@graphs), collapse = ", "))
graph_name <- snn_graphs[1]
message(sprintf("  Using graph: %s", graph_name))
so <- FindClusters(so, graph.name = graph_name, resolution = RESOLUTION, verbose = FALSE)
message(sprintf("  Clusters: %d", length(unique(so$seurat_clusters))))
gc()


################################################################################
# STEP 5: RunUMAP on X_scVI
################################################################################

message("RunUMAP on X_scVI embedding ...")
so <- RunUMAP(so, reduction = "X_scVI", dims = SCVI_DIMS,
              reduction.name = "umap_scvi", verbose = FALSE)
message("  UMAP done")


################################################################################
# STEP 6: Visualise
################################################################################

message("Generating UMAP plots ...")

make_plot <- function(grp, title, label = FALSE) {
  if (!grp %in% colnames(so@meta.data)) return(NULL)
  DimPlot(so, group.by = grp, reduction = "umap_scvi",
          label = label, repel = label, raster = TRUE) + ggtitle(title)
}

plots <- Filter(Negate(is.null), list(
  make_plot("seurat_clusters", "Clusters",      label = TRUE),
  make_plot("condition",       "By Condition"),
  make_plot("sample_id",       "By Sample"),
  make_plot("technology",      "By Technology"),
  make_plot("leiden_scVI",     "scVI Leiden")
))

n_cols <- min(2L, length(plots))
pdf(file.path(OUT_DIR, "RenalUrothelium_UMAP_overview.pdf"),
    width = 18, height = ceiling(length(plots) / n_cols) * 7)
print(wrap_plots(plots, ncol = n_cols))
dev.off()
message("  Saved: RenalUrothelium_UMAP_overview.pdf")

for (ct_col in c("cell_type_scanvi", "cell_type_original")) {
  if (ct_col %in% colnames(so@meta.data)) {
    so@meta.data[[ct_col]] <- as.character(so@meta.data[[ct_col]])
    so@meta.data[[ct_col]][is.na(so@meta.data[[ct_col]])] <- "Unknown"
    p_ct <- DimPlot(
      so, group.by = ct_col, reduction = "umap_scvi",
      label = TRUE, repel = TRUE, raster = TRUE
    ) + ggtitle(ct_col)
    pdf(file.path(OUT_DIR, paste0("RenalUrothelium_", ct_col, ".pdf")),
        width = 14, height = 10)
    print(p_ct)
    dev.off()
    message(sprintf("  Saved: RenalUrothelium_%s.pdf", ct_col))
  }
}

markers <- c("Krt5", "Krt14", "Krt20", "Krt8", "Krt18",
             "Upk2", "Upk3a", "Trp63", "Upk1a", "Upk1b")
present <- intersect(markers, rownames(so))

if (length(present) > 0) {
  pdf(file.path(OUT_DIR, "RenalUrothelium_uro_markers.pdf"),
      width = 18, height = 8)
  print(FeaturePlot(so, features = present, reduction = "umap_scvi",
                    ncol = 5, raster = TRUE))
  dev.off()
  message("  Saved: RenalUrothelium_uro_markers.pdf")
}


################################################################################
# STEP 7: Save
################################################################################

message(sprintf("Saving to %s ...", OUT_PATH))
saveRDS(so, OUT_PATH)

message("\n===== Renal Urothelium scVI clustering complete =====")
message(sprintf("Cells    : %s", format(ncol(so), big.mark = ",")))
message(sprintf("Clusters : %d", length(unique(so$seurat_clusters))))
message(sprintf("Output   : %s", OUT_PATH))

