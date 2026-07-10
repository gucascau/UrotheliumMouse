################################################################################
# 03_Xenium_Deconvolution.R
#
# RCTD cell-type deconvolution of 12 Xenium UUO samples using the annotated
# kidney scRNA-seq as reference (final_annotation_with_uro).
#
# Strategy mirrors 02_VisiumHD_Deconvolution.r:
#   1. Build RCTD reference from annotated kidney scRNA-seq
#   2. Loop over each sample -> build SpatialRNA query from that sample's FOV
#   3. Run RCTD per sample (doublet mode: Xenium cells are single-cell resolution)
#   4. Combine weights across all samples, add to Seurat object
#   5. Save and plot
#
# Input:  Xenium_harmony_integrated_prelabeltransfer.rds
#         RenalUrothelium_allcells_scvi_annotations_meta_uro_updated.rds
# Output: Xenium_RCTD_deconvolved.rds
#         Xenium_RCTD_dominant_celltype.pdf
################################################################################

library(Seurat)
library(spacexr)
library(ggplot2)
library(dplyr)

options(future.globals.maxSize = 16 * 1024^3)

log_mem <- function(label) {
  g  <- gc(verbose = FALSE)
  mb <- sum(g[, 2])
  message(sprintf("  [mem] %s: %.1f GB in use", label, mb / 1024))
}

# ── Paths ─────────────────────────────────────────────────────────────────────
OUT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/SpatialScripts/output"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OBJECT_RDS <- file.path(OUT_DIR, "Xenium_harmony_integrated_prelabeltransfer.rds")
REF_RDS    <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/RenalUrotheliumScripts/output/RenalUrothelium_allcells_scvi_annotations_meta_uro_updated.rds"
OUT_DECONV <- file.path(OUT_DIR, "Xenium_RCTD_deconvolved.rds")
OUT_DECONV_UPdated <- file.path(OUT_DIR, "Xenium_RCTD_updated_deconvolved.rds")
OUT_PDF    <- file.path(OUT_DIR, "Xenium_RCTD_dominant_celltype.pdf")

# Sample IDs (must match sample_id metadata column in object)
SAMPLE_IDS <- c("ShamL","ShamR","Hour4L","Hour4R",
                 "Hour12L","Hour12R","Day2L","Day2R",
                 "Day14L","Day14R","Week6L","Week6R")

# # ── Load Xenium integrated object ─────────────────────────────────────────────
# message("==> Loading Xenium integrated object ...")
# object <- readRDS(OBJECT_RDS)
# object <- JoinLayers(object, assay = "Xenium")
# log_mem("after loading object")

# # Report available FOV slots (used for per-sample coordinate lookup)
# all_fovs <- Images(object)
# message("Available FOV slots: ", paste(all_fovs, collapse = ", "))

# # ── Build RCTD reference from annotated kidney scRNA-seq ──────────────────────
# message("==> Building RCTD reference ...")
# ref_so <- readRDS(REF_RDS)
# log_mem("after loading ref_so")

# ref_counts <- GetAssayData(ref_so, assay = "RNA", layer = "counts")
# cell_types  <- droplevels(as.factor(ref_so$final_annotation_with_uro))
# names(cell_types) <- colnames(ref_so)
# nUMI <- colSums(ref_counts)

# rctd_ref <- spacexr::Reference(
#   counts     = ref_counts,
#   cell_types = cell_types,
#   nUMI       = nUMI,
#   min_UMI    = 1,
#   n_max_cells  = Inf
# )
# rm(ref_so, ref_counts, nUMI, cell_types); gc()
# log_mem("after building RCTD reference")

# message("Reference dimensions (genes x cells): ", paste(dim(rctd_ref@counts), collapse = " x "))
# message("Cell type distribution:")
# print(table(rctd_ref@cell_types))

# # ── Run RCTD per sample ───────────────────────────────────────────────────────
# # Each Xenium sample has its own coordinate space in its own FOV slot.
# # Running RCTD per-sample avoids coordinate collisions across samples.
# # doublet mode: Xenium provides single-cell resolution — each cell is at most
# # a mixture of 2 types (imperfect segmentation boundary cells).
# all_rctd_meta <- list()

# for (i in seq_along(SAMPLE_IDS)) {
#   samp <- SAMPLE_IDS[i]
#   message(sprintf("\n==> RCTD for sample: %s", samp))

#   # FOV slots are all named "fov" at LoadXenium() time (one per sample); merge()
#   # renumbers the collisions as fov, fov.2, fov.3, ... in merge order, which is
#   # the same order as SAMPLE_IDS. Verified by cross-tabulating each FOV's cell
#   # barcodes against the sample_id metadata column (see
#   # verify_fov_sample_mapping.R): each FOV slot maps 1:1 to exactly one sample,
#   # in SAMPLE_IDS order.
#   fov_name <- all_fovs[i]

#   # Cells belonging to this sample
#   samp_cells <- WhichCells(object, expression = sample_id == samp)
#   message(sprintf("  Cells: %d", length(samp_cells)))

#   # Spatial coordinates from this sample's FOV
#   fov_coords <- GetTissueCoordinates(object, image = fov_name, which = "centroids")
#   # fov_coords columns: cell, x, y — keep only cells present in samp_cells
#   fov_coords <- fov_coords[fov_coords$cell %in% samp_cells, ]

#   coords_df <- data.frame(
#     x = fov_coords$x,
#     y = fov_coords$y,
#     row.names = fov_coords$cell
#   )

#   # Count matrix subset to this sample
#   samp_counts <- GetAssayData(object, assay = "Xenium", layer = "counts")[, fov_coords$cell]

#   # Gene overlap: Xenium gene names may have dashes instead of underscores
#   xenium_genes <- rownames(samp_counts)
#   ref_genes    <- rownames(rctd_ref@counts)
#   shared_genes <- intersect(xenium_genes, ref_genes)
#   message(sprintf("  Shared genes with reference: %d / %d Xenium genes",
#                   length(shared_genes), length(xenium_genes)))

#   samp_counts <- samp_counts[shared_genes, ]
#   samp_nUMI   <- colSums(samp_counts)

#   query <- SpatialRNA(coords_df, samp_counts, samp_nUMI)
#   rm(fov_coords, coords_df, samp_counts, samp_nUMI); gc()

#   # Run RCTD
#   myRCTD <- create.RCTD(query, rctd_ref, max_cores = 20, CELL_MIN_INSTANCE = 1)
#   rm(query); gc()
#   myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")

#   # Extract and store weights
#   weights <- as.matrix(normalize_weights(myRCTD@results$weights))
#   colnames(weights) <- paste0("rctd_", colnames(weights))

#   meta_i <- as.data.frame(weights)
#   meta_i$rctd_dominant_celltype <- sub(
#     "^rctd_", "",
#     colnames(weights)[max.col(weights, ties.method = "first")]
#   )

#   all_rctd_meta[[samp]] <- meta_i
#   rm(myRCTD, weights, meta_i); gc()
#   log_mem(paste("after RCTD for", samp))
# }
# # save the all_rectd_meta list for debugging
# saveRDS(all_rctd_meta, file.path(OUT_DIR, "all_rctd_meta.rds"))

# # ── Combine weights across all samples and add to Seurat object ───────────────
# message("\n==> Adding RCTD weights to Seurat object ...")
# combined_meta <- do.call(rbind, unname(all_rctd_meta))
# rm(all_rctd_meta, rctd_ref); gc()

# # AddMetaData matches on rownames — cells not in combined_meta get NA
# object <- AddMetaData(object, combined_meta)
# rm(combined_meta); gc()
# log_mem("after adding metadata")

# message("Dominant cell type distribution:")
# print(table(object$rctd_dominant_celltype, useNA = "ifany"))

# # ── Save ──────────────────────────────────────────────────────────────────────
# message("==> Saving deconvolved object ...")
# saveRDS(object, OUT_DECONV)
# message("  Saved: ", OUT_DECONV)

# # ── UMAP colored by dominant cell type ────────────────────────────────────────
# message("==> Plotting ...")
# pdf(OUT_PDF, width = 14, height = 6)

# p_umap <- DimPlot(object,
#   reduction = "umap", group.by = "rctd_dominant_celltype",
#   label = TRUE, repel = TRUE, pt.size = 0.1
# ) + ggtitle("RCTD dominant cell type (UMAP)")

# p_clust <- DimPlot(object,
#   reduction = "umap", group.by = "seurat_clusters",
#   label = TRUE, repel = TRUE, pt.size = 0.1
# ) + ggtitle("Seurat clusters (UMAP)")

# print(p_umap | p_clust)

# # Spatial plot per sample
# for (i in seq_along(SAMPLE_IDS)) {
#   samp     <- SAMPLE_IDS[i]
#   fov_name <- Images(object)[i]

#   samp_cells <- WhichCells(object, expression = sample_id == samp)
#   sub_obj    <- subset(object, cells = samp_cells)

#   tryCatch({
#     print(
#       ImageDimPlot(sub_obj, fov = fov_name, group.by = "rctd_dominant_celltype",
#                    cols = "polychrome", axes = TRUE, size = 0.5) +
#         ggtitle(paste("RCTD cell types:", samp))
#     )
#   }, error = function(e) {
#     message("  ImageDimPlot failed for ", samp, ": ", e$message)
#   })
# }

# dev.off()
# message("  Saved: ", OUT_PDF)

# message("==> Done.")
# log_mem("final")


# # Due to 300 genes in the Xenium Panel, we cannot fully capture the cell types for all these cells, therefore we count on the biomarker to manually annoate these cells

object<- readRDS(OUT_DECONV)


# the previous clusters are not well categorized into different cell types, I further increased the resolution. 
object <- FindClusters(object, resolution = 1, cluster.name = "seurat_clusters")
object <- RunUMAP(object,
  reduction      = "harmony",
  dims           = 1:30,
  reduction.name = "umap"
)

# save the reclustered object
saveRDS(object, OUT_DECONV_UPdated)

object<- readRDS(OUT_DECONV_UPdated)
# check the Identity of the object
Idents(object)


markers <- FindAllMarkers(
  object,
  only.pos = TRUE,      # only upregulated markers, not both directions
  min.pct = 0.25,        # gene must be detected in ≥25% of cells in at least one group
  logfc.threshold = 0.25 # minimum log2FC to bother testing
)

markers %>% pull(cluster) %>% table()


library(dplyr)
# save the biomarkers
saveRDS(markers, file.path(OUT_DIR, "markers.rds"))
markers<- readRDS(file.path(OUT_DIR, "markers.rds"))

top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10)


write.csv(top_markers, file.path(OUT_DIR, "Xenium_clusters_topmarker.csv"))
TopMarkersHeatmap<- DoHeatmap(object , features = top_markers$gene)
# save the heatmap
ggsave("Xenium_TopMarkersHeatmap.pdf", plot= TopMarkersHeatmap, height = 10, width = 10)

# ── Manual cluster annotation ──────────────────────────────────────────────
# CSV has a UTF-8 BOM that can mangle read.csv's column names, so match by
# position: col 1 = Cluster, col 2 = Cell Type, col 3 = Abbreviation,
# col 4 = Representative Markers.
ClusterType <- read.csv(file.path(OUT_DIR, "Xenium_clusterannotation.csv"))

cluster_to_celltype <- setNames(ClusterType[[2]], as.character(ClusterType[[1]]))
cluster_to_abbr     <- setNames(ClusterType[[3]], as.character(ClusterType[[1]]))
cluster_to_repgenes <- setNames(ClusterType[[4]], as.character(ClusterType[[1]]))
cluster_to_manconf <- setNames(ClusterType[[5]], as.character(ClusterType[[1]]))


object$CellType         <- unname(cluster_to_celltype[as.character(object$seurat_clusters)])
object$CellTypeAbbr     <- unname(cluster_to_abbr[as.character(object$seurat_clusters)])
object$CellTypeRepGenes <- unname(cluster_to_repgenes[as.character(object$seurat_clusters)])
object$CellTypeManualAnnConfidence <- unname(cluster_to_manconf[as.character(object$seurat_clusters)])
#

# ── Reconcile manual annotation with RCTD per-cell calls ───────────────────
# Manual review recorded a per-cluster confidence in the CSV
# (CellTypeManualAnnConfidence: High/Low):
#   - High confidence -> always use the manual CellTypeAbbr.
#   - Low confidence  -> use RCTD's own per-cell call instead, since manual
#     review flagged the cluster as ambiguous/mixed; if RCTD has no call for
#     that cell (NA), fall back to the manual (low-confidence) label, tagged
#     with an "_LQ" suffix so these low-quality fallback calls stay
#     distinguishable from a confident manual annotation.
# A few RCTD labels never override the manual label even in low-confidence
# clusters — vascular subtypes the manual clustering already groups under
# Peri/VSMC or Endo. (Note: the RCTD reference spells it "Vas-Afferens",
# not "Vas-Afference" — fixed here.)
rctd_rescue_exclude <- c("Per", "Vas-Efferens", "Vas-Afferens",
                          "Glom-Endo", "Desc-Vasa-Recta", "Asc-Vasa-Recta")

# RCTD reference labels that name the same cell type as a manual abbreviation
# but with different spelling/style get recoded to the manual spelling, so
# FinalCellType never has two strings for one cell type. CTAL/MTAL/LOH are
# all folded into "TAL" per user confirmation; PTS1/PTS3 map 1:1 onto the
# existing PT-S1/PT-S3 abbreviations (PTS2/PTS3T2 are left unrenamed since
# folding them into the combined "PT-S2/S3" bucket is a many-to-one call the
# user chose not to make).
rctd_name_map <- c(
  "Urothelium" = "Uro",
  "Macro"      = "Mac",
  "PTS1"       = "PT-S1",
  "PTS3"       = "PT-S3",
  "ICA"        = "A-IC",
  "ICB"        = "B-IC",
  "Per"        = "Peri/VSMC",
  "CTAL"       = "TAL",
  "MTAL"       = "TAL",
  "LOH"        = "TAL"
)

is_low_conf <- object$CellTypeManualAnnConfidence == "Low"
rctd_call   <- object$rctd_dominant_celltype
recode_idx  <- rctd_call %in% names(rctd_name_map)
rctd_call[recode_idx] <- rctd_name_map[rctd_call[recode_idx]]

use_rctd_idx <- is_low_conf & !is.na(rctd_call) & !(rctd_call %in% rctd_rescue_exclude)
na_fallback_idx <- is_low_conf & is.na(rctd_call)

final_celltype <- object$CellTypeAbbr
final_celltype[use_rctd_idx] <- rctd_call[use_rctd_idx]
object$FinalCellType <- final_celltype

# Flag cells that fell back to the manual (low-confidence) label because RCTD
# had no call for them, in a separate column rather than mangling the label
# itself with a suffix.
object$FinalCellTypeFlag <- NA_character_
object$FinalCellTypeFlag[na_fallback_idx] <- "LQ"

message("Final cell type distribution:")
print(table(object$FinalCellType, useNA = "ifany"))
message("Low-quality (RCTD-NA) fallback cells: ", sum(na_fallback_idx))

# ── Save annotated object ───────────────────────────────────────────────────
saveRDS(object, OUT_DECONV_UPdated)
message("Saved annotated object: ", OUT_DECONV_UPdated)
