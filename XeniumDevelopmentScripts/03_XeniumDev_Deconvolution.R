################################################################################
# 03_XeniumDev_Deconvolution.R
#
# RCTD cell-type deconvolution of the 4 Xenium developmental-kidney samples
# (GSE286051: W12/W92, male/female) using the adult RenalUrothelium reference
# (final_annotation_with_uro) -- same reference and "doublet" mode as
# SpatialScripts/03_Xenium_Updated_Deconvolution.R's 12-sample UUO pipeline
# (Xenium is near-single-cell resolution, so doublet mode -- at most 2 cell
# types per cell from imperfect segmentation -- not "full" mode).
#
# Unlike VisiumLowDevelopmentScripts' Chen2025 stage-matched deconvolution,
# these 4 samples are adult (W12, W92) -- the adult RenalUrothelium reference
# is appropriate here without age-matching concerns.
#
# IMPORTANT (memory from VisiumLowDiseaseScripts/04_VisiumLowDisease_
# Deconvolution.R): building an RCTD reference from this same 968,373-cell
# RenalUrothelium object and calling create.RCTD() OOM'd at 200G/16 cores --
# the established precedent for this exact reference
# (SpatialScripts/submit_02_VisiumHD_deconvolution.sh,
# submit_03_Xenium_deconvolution.sh) requests 800-900G+ and 24-30 cores.
# submit_03_XeniumDev_Deconvolution.sh is sized to match that precedent from
# the start.
#
# spacexr::Reference() rejects "/" in cell type names (hit this with
# Chen2025's "B/Plasma" label) -- sanitized defensively even though
# final_annotation_with_uro doesn't currently contain one.
#
# rctd_dominant_celltype requires the winning RCTD weight to be a true
# MAJORITY (>0.5 of that cell's total weight), not just a plurality --
# max.col() alone would call a cell type "dominant" even at, say, 30%
# weight if nothing else scored higher. Confirmed interactively this
# matters in practice for Urothelium specifically (median winning weight
# was 0.59, 36.5% below 0.5, dominant contaminant Asc-Vasa-Recta present in
# 37% of "Urothelium" calls -- a real adjacent-anatomy confound, not noise).
# Cells that don't clear the bar get rctd_dominant_celltype = NA rather
# than an unreliable label.
#
# Input:  XeniumDevelopmentScripts/output/XeniumDev_harmony_integrated_prelabeltransfer.rds
#         RenalUrotheliumScripts/output/RenalUrothelium_allcells_scvi_annotations_meta_uro_updated.rds
# Output: XeniumDev_RCTD_deconvolved.rds
#         XeniumDev_RCTD_dominant_celltype.pdf
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(spacexr)
  library(ggplot2)
  library(dplyr)
})

options(future.globals.maxSize = 32 * 1024^3)

log_mem <- function(label) {
  g  <- gc(verbose = FALSE)
  mb <- sum(g[, 2])
  message(sprintf("  [mem] %s: %.1f GB in use", label, mb / 1024))
}

SCRIPT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/XeniumDevelopmentScripts"
OUT_DIR    <- file.path(SCRIPT_DIR, "output")
OBJECT_RDS <- file.path(OUT_DIR, "XeniumDev_harmony_integrated_prelabeltransfer.rds")
REF_RDS    <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/RenalUrotheliumScripts/output/RenalUrothelium_allcells_scvi_annotations_meta_uro_updated.rds"
OUT_DECONV <- file.path(OUT_DIR, "XeniumDev_RCTD_deconvolved.rds")
OUT_PDF    <- file.path(OUT_DIR, "XeniumDev_RCTD_dominant_celltype.pdf")

# ── Load Xenium integrated object ───────────────────────────────────────────
message("==> Loading Xenium integrated object ...")
if (!file.exists(OBJECT_RDS)) {
  stop("Missing ", OBJECT_RDS, " -- run 02_XeniumDev_Integrate_Harmony.R first.")
}
object <- readRDS(OBJECT_RDS)
object <- JoinLayers(object, assay = "Xenium")
log_mem("after loading object")

all_fovs <- Images(object)
message("Available FOV slots: ", paste(all_fovs, collapse = ", "))

# ── Build RCTD reference from annotated kidney scRNA-seq ────────────────────
message("==> Building RCTD reference ...")
ref_so <- readRDS(REF_RDS)
log_mem("after loading ref_so")

ref_counts <- GetAssayData(ref_so, assay = "RNA", layer = "counts")
# spacexr::Reference() rejects "/" in cell type names -- see header comment.
cell_types_chr <- gsub("[/]", "_", as.character(ref_so$final_annotation_with_uro))
cell_types <- droplevels(as.factor(cell_types_chr))
names(cell_types) <- colnames(ref_so)
nUMI <- colSums(ref_counts)

rctd_ref <- spacexr::Reference(
  counts     = ref_counts,
  cell_types = cell_types,
  nUMI       = nUMI,
  min_UMI    = 1,
  n_max_cells  = Inf
)
rm(ref_so, ref_counts, nUMI, cell_types); gc()
log_mem("after building RCTD reference")

message("Reference dimensions (genes x cells): ", paste(dim(rctd_ref@counts), collapse = " x "))
message("Cell type distribution:")
print(table(rctd_ref@cell_types))

# ── Run RCTD per sample (doublet mode -- Xenium is single-cell resolution) ──
all_rctd_meta <- list()

for (i in seq_along(all_fovs)) {
  fov_name <- all_fovs[i]
  fov_coords <- GetTissueCoordinates(object, image = fov_name, which = "centroids")
  samp_cells <- fov_coords$cell
  samp <- unique(object$sample_id[samp_cells])
  stopifnot(length(samp) == 1)
  message(sprintf("\n==> RCTD for sample: %s (fov %s, %d cells)", samp, fov_name, length(samp_cells)))

  coords_df <- data.frame(x = fov_coords$x, y = fov_coords$y, row.names = fov_coords$cell)

  samp_counts <- GetAssayData(object, assay = "Xenium", layer = "counts")[, samp_cells]
  shared_genes <- intersect(rownames(samp_counts), rownames(rctd_ref@counts))
  message(sprintf("  Shared genes with reference: %d / %d Xenium genes",
                  length(shared_genes), nrow(samp_counts)))
  samp_counts <- samp_counts[shared_genes, ]
  samp_nUMI <- colSums(samp_counts)

  query <- SpatialRNA(coords_df, samp_counts, samp_nUMI)
  rm(fov_coords, coords_df, samp_counts, samp_nUMI); gc()

  myRCTD <- create.RCTD(query, rctd_ref, max_cores = 16, CELL_MIN_INSTANCE = 1)
  rm(query); gc()
  myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")

  weights <- as.matrix(normalize_weights(myRCTD@results$weights))
  colnames(weights) <- paste0("rctd_", colnames(weights))

  meta_i <- as.data.frame(weights)
  # normalize_weights() makes each cell's weights sum to 1 across at most 2
  # cell types (doublet mode), so max.col() alone only guarantees a
  # PLURALITY winner, not a majority -- confirmed interactively on the
  # Urothelium call specifically that this matters: median winning weight
  # was only 0.59, and 36.5% of "Urothelium" cells were below 0.5, with the
  # dominant contaminant (Asc-Vasa-Recta, a real adjacent-anatomy confound)
  # present in 37% of them. Requiring MIN_DOMINANT_WEIGHT here applies that
  # same majority-not-plurality bar to all 36 cell types, not just
  # Urothelium -- cells whose top RCTD weight doesn't clear it get NA
  # rather than a low-confidence label.
  MIN_DOMINANT_WEIGHT <- 0.5
  top_weight <- apply(weights, 1, max)
  meta_i$rctd_dominant_celltype <- ifelse(
    top_weight > MIN_DOMINANT_WEIGHT,
    sub("^rctd_", "", colnames(weights)[max.col(weights, ties.method = "first")]),
    NA_character_
  )

  all_rctd_meta[[samp]] <- meta_i
  rm(myRCTD, weights, meta_i); gc()
  log_mem(paste("after RCTD for", samp))
}

# ── Combine weights across all samples and add to Seurat object ────────────
message("\n==> Adding RCTD weights to Seurat object ...")
combined_meta <- do.call(rbind, unname(all_rctd_meta))
rm(all_rctd_meta, rctd_ref); gc()

object <- AddMetaData(object, combined_meta)
rm(combined_meta); gc()
log_mem("after adding metadata")

message("Dominant cell type distribution:")
print(table(object$rctd_dominant_celltype, useNA = "ifany"))

# ── Save ─────────────────────────────────────────────────────────────────────
message("==> Saving deconvolved object ...")
saveRDS(object, OUT_DECONV)
message("  Saved: ", OUT_DECONV)

# ── UMAP + per-sample spatial plot of dominant cell type ────────────────────
message("==> Plotting ...")
pdf(OUT_PDF, width = 14, height = 6)

p_umap <- DimPlot(object, reduction = "umap", group.by = "rctd_dominant_celltype",
  label = TRUE, repel = TRUE, pt.size = 0.1) + ggtitle("RCTD dominant cell type (UMAP)")
p_clust <- DimPlot(object, reduction = "umap", group.by = "seurat_clusters",
  label = TRUE, repel = TRUE, pt.size = 0.1) + ggtitle("Seurat clusters (UMAP)")
print(p_umap | p_clust)

for (i in seq_along(all_fovs)) {
  fov_name <- all_fovs[i]
  fov_coords <- GetTissueCoordinates(object, image = fov_name, which = "centroids")
  samp <- unique(object$sample_id[fov_coords$cell])
  sub_obj <- subset(object, cells = fov_coords$cell)
  tryCatch({
    print(
      ImageDimPlot(sub_obj, fov = fov_name, group.by = "rctd_dominant_celltype",
                   cols = "polychrome", axes = TRUE, size = 0.5) +
        ggtitle(paste("RCTD cell types:", samp))
    )
  }, error = function(e) {
    message("  ImageDimPlot failed for ", samp, ": ", e$message)
  })
}

dev.off()
message("  Saved: ", OUT_PDF)

message("==> Done.")
log_mem("final")
