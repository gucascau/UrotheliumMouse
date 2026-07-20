################################################################################
# 04_XeniumDev_MarkerHeatmap.R
#
# Cluster marker genes + a blank annotation-template CSV for the 4 Xenium
# developmental-kidney samples, mirroring the marker/heatmap half of
# SpatialScripts/03_Xenium_Updated_Deconvolution.R -- but stopping there.
# That script's actual cluster->cell-type calls (Xenium_clusterannotation.csv:
# CellType/Abbreviation/RepresentativeMarkers/ManualConfidence per cluster)
# came from a human reviewing marker genes against known kidney biology,
# which this script cannot substitute for. It writes the same shape of CSV
# template (with the FindAllMarkers top genes pre-filled for reference) with
# the annotation columns left blank for manual review, rather than guessing.
#
# Input:  XeniumDevelopmentScripts/output/XeniumDev_RCTD_deconvolved.rds
# Output: XeniumDev_clusters_topmarkers.csv
#         XeniumDev_TopMarkersHeatmap.pdf
#         XeniumDev_clusterannotation_template.csv  (fill in by hand, then
#           feed into whatever script reconciles RCTD + manual calls, mirroring
#           the "Manual cluster annotation" section of 03_Xenium_Updated_
#           Deconvolution.R)
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

SCRIPT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/XeniumDevelopmentScripts"
OUT_DIR    <- file.path(SCRIPT_DIR, "output")
OBJECT_RDS <- file.path(OUT_DIR, "XeniumDev_RCTD_deconvolved.rds")

message("==> Loading deconvolved Xenium object ...")
if (!file.exists(OBJECT_RDS)) {
  stop("Missing ", OBJECT_RDS, " -- run 03_XeniumDev_Deconvolution.R first.")
}
object <- readRDS(OBJECT_RDS)

message("==> Finding cluster marker genes ...")
markers <- FindAllMarkers(
  object,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25
)
saveRDS(markers, file.path(OUT_DIR, "XeniumDev_markers.rds"))

message("Cells per cluster:")
print(table(object$seurat_clusters))

top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10)
write.csv(top_markers, file.path(OUT_DIR, "XeniumDev_clusters_topmarkers.csv"), row.names = FALSE)
message("  Saved: ", file.path(OUT_DIR, "XeniumDev_clusters_topmarkers.csv"))

message("==> Plotting top-marker heatmap ...")
p_heatmap <- DoHeatmap(object, features = top_markers$gene)
ggplot2::ggsave(file.path(OUT_DIR, "XeniumDev_TopMarkersHeatmap.pdf"), plot = p_heatmap, height = 10, width = 10)
message("  Saved: ", file.path(OUT_DIR, "XeniumDev_TopMarkersHeatmap.pdf"))

# ── Blank annotation template (for manual review) ───────────────────────────
# One row per cluster, pre-filled with its top marker genes and RCTD's most
# common per-cell call in that cluster (a starting hint, not a final answer)
# -- CellType/Abbreviation/ManualConfidence left blank for hand annotation.
message("==> Writing cluster annotation template ...")
clusters <- sort(unique(object$seurat_clusters))
rctd_mode_by_cluster <- sapply(clusters, function(cl) {
  calls <- object$rctd_dominant_celltype[object$seurat_clusters == cl]
  calls <- calls[!is.na(calls)]
  if (length(calls) == 0) return(NA_character_)
  names(sort(table(calls), decreasing = TRUE))[1]
})
top_genes_by_cluster <- sapply(clusters, function(cl) {
  paste(top_markers$gene[top_markers$cluster == cl], collapse = ", ")
})

template <- data.frame(
  Cluster                    = clusters,
  CellType                   = "",
  Abbreviation               = "",
  RepresentativeMarkers      = top_genes_by_cluster,
  ManualConfidence           = "",
  RCTD_mode_celltype_hint    = rctd_mode_by_cluster,
  stringsAsFactors = FALSE
)
out_template <- file.path(OUT_DIR, "XeniumDev_clusterannotation_template.csv")
write.csv(template, out_template, row.names = FALSE)
message("  Saved: ", out_template)
message("\n  Fill in CellType/Abbreviation/ManualConfidence by hand (see marker")
message("  heatmap + RCTD hint column), then reconcile with RCTD per-cell calls")
message("  the way 03_Xenium_Updated_Deconvolution.R's manual-annotation section")
message("  does, if that's the analysis you want next.")

message("\n==> Done.")
