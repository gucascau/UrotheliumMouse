################################################################################
# 01c_VisiumHD_full_FindAllMarkers.R
#
# FindAllMarkers on the finished bladder + 2-kidney Harmony/sketch integration
# (seurat_cluster.harmony.projected, 15 clusters), mirroring the FindAllMarkers
# call used for the kidney-only object in 01_VisiumHD_integrate_harmony_sketch.R.
#
# Input : output/VisiumHD_harmony_integrated_clustered.rds
# Output: output/UnBiasMarkersSpatialHD_Full.rds
################################################################################

.libPaths(c("/home/gdbecknelllab/xxw004/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(MAST)
})

options(future.globals.maxSize = 16 * 1024^3)

OUT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/SpatialScripts/output"
IN_RDS  <- file.path(OUT_DIR, "VisiumHD_harmony_integrated_clustered.rds")
OUT_RDS <- file.path(OUT_DIR, "UnBiasMarkersSpatialHD_Full.rds")

message("==> Loading: ", IN_RDS)
object <- readRDS(IN_RDS)
DefaultAssay(object) <- "Spatial.008um"
Idents(object) <- "seurat_cluster.harmony.projected"

message("==> JoinLayers ...")
object <- JoinLayers(object)

message("==> FindAllMarkers (MAST) ...")
UnBiasMarkersSpatialHD_Full <- FindAllMarkers(object, test.use = "MAST", verbose = TRUE)

saveRDS(UnBiasMarkersSpatialHD_Full, file = OUT_RDS)
message("  Saved: ", OUT_RDS)
message("==> Done.")
