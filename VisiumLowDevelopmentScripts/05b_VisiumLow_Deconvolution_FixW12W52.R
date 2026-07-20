################################################################################
# 05b_VisiumLow_Deconvolution_FixW12W52.R
#
# Targeted fix-up: 05_VisiumLow_Deconvolution.R's original run left all 4
# W12 and W52 samples with NA in every rctd_* column -- create.RCTD() failed
# for these two stages only ("'dims' cannot be of length 0", from RCTD's
# get_cell_mean() choking on the W12/W52 Chen2025 references' singleton/
# near-singleton cell types: LOH_AL_proliferating n=1, UBP n=3), silently
# swallowed by that script's per-sample tryCatch. Root cause and fix (drop
# cell types with < MIN_CELLS_PER_TYPE reference cells before building
# spacexr::Reference(), not just via create.RCTD()'s CELL_MIN_INSTANCE) are
# now folded into 05_VisiumLow_Deconvolution.R itself for any future full
# re-run; this script re-runs RCTD for only the 4 affected samples against
# the corrected references and patches VisiumLow_deconvolved.rds in place,
# rather than repeating the full 24-sample, multi-hour run.
#
# Affected samples:
#   W12: GSM8704853 (AKICT3F), GSM8704854 (AKICT3M)
#   W52: GSM8704844 (NMK52-F), GSM8704845 (NMK52-M)
#
# Input:  VisiumLowDevelopmentScripts/output/VisiumLow_deconvolved.rds
#         (written by 05_VisiumLow_Deconvolution.R)
#         UsedSpatialData/VisiumLow/DevelopedKidney/*_obj.rds (raw per-sample
#         objects, re-loaded fresh rather than pulling from the merged
#         object, matching 05's own loading logic)
#         Chen2025 raw counts + metadata (same sources as 05)
# Output: VisiumLow_deconvolved.rds (patched in place, rctd_* columns filled
#           in for the 4 previously-NA samples)
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(spacexr)
  library(dplyr)
})

options(future.globals.maxSize = 16 * 1024^3)

log_mem <- function(label) {
  g  <- gc(verbose = FALSE)
  mb <- sum(g[, 2])
  message(sprintf("  [mem] %s: %.1f GB in use", label, mb / 1024))
}

DATA_DIR   <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSpatialData/VisiumLow/DevelopedKidney"
SCRIPT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/VisiumLowDevelopmentScripts"
OUT_DIR    <- file.path(SCRIPT_DIR, "output")
CHEN_RDS   <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/RawMouseSingleCellDatasets/MultiOmicSpatialMouseKidney_SexSpecific_lifespan_Chen2025NatGenet_zellkonvertedConverted.rds"
CHEN_COUNTS_DIR <- file.path(OUT_DIR, "Chen2025_rawcounts")
OBJECT_RDS <- file.path(OUT_DIR, "VisiumLow_deconvolved.rds")

MIN_CELLS_PER_TYPE <- 25  # see 05_VisiumLow_Deconvolution.R for rationale

FIX_SAMPLES <- data.frame(
  rds_file = c("GSM8704853_20210129-AKICT3F-Fc1U1Z1Bs1_obj.rds",
               "GSM8704854_20210129-AKICT3M-Fc1U1Z1Bs1_obj.rds",
               "GSM8704844_20201207-NMK52-F-Fc1U1Z1Bs1_obj.rds",
               "GSM8704845_20201207-NMK52-M-Fc1U1Z1Bs1_obj.rds"),
  sample_id = c("GSM8704853_20210129-AKICT3F-Fc1U1Z1Bs1",
                "GSM8704854_20210129-AKICT3M-Fc1U1Z1Bs1",
                "GSM8704844_20201207-NMK52-F-Fc1U1Z1Bs1",
                "GSM8704845_20201207-NMK52-M-Fc1U1Z1Bs1"),
  stage = c("W12", "W12", "W52", "W52"),
  stringsAsFactors = FALSE
)

# ── Load Chen2025 reference (metadata + raw counts, same as 05) ────────────
message("==> Loading Chen2025 reference metadata + raw counts ...")
chen_meta_obj <- readRDS(CHEN_RDS)
chen_meta <- chen_meta_obj@meta.data[, c("Age", "celltype_final")]
rm(chen_meta_obj); gc()

chen_counts <- Read10X(data.dir = CHEN_COUNTS_DIR, unique.features = TRUE)
chen_meta$Age <- factor(chen_meta$Age, levels = c("E16.5", "P0", "W3", "W12", "W52", "W92"))
chen_meta$celltype_final <- gsub("[/ ]", "_", as.character(chen_meta$celltype_final))
chen_meta$celltype_final <- droplevels(as.factor(chen_meta$celltype_final))
ALL_CELLTYPES <- levels(chen_meta$celltype_final)
log_mem("after loading Chen2025 reference")

# ── Build corrected references for just W12 and W52 ────────────────────────
message("\n==> Building corrected W12 and W52 RCTD references ...")
stage_refs <- list()
for (stage in c("W12", "W52")) {
  stage_cells_all <- rownames(chen_meta)[chen_meta$Age == stage]
  stage_types_all <- droplevels(chen_meta[stage_cells_all, "celltype_final"])
  type_counts <- table(stage_types_all)
  dropped_types <- names(type_counts)[type_counts < MIN_CELLS_PER_TYPE]
  message(sprintf("  %s: dropping cell type(s) with < %d reference cells: %s",
                   stage, MIN_CELLS_PER_TYPE, paste(dropped_types, collapse = ", ")))
  stage_cells <- stage_cells_all[!stage_types_all %in% dropped_types]

  stage_counts <- chen_counts[, stage_cells]
  stage_types  <- droplevels(chen_meta[stage_cells, "celltype_final"])
  names(stage_types) <- stage_cells
  stage_nUMI <- colSums(stage_counts)

  stage_refs[[stage]] <- spacexr::Reference(
    counts = stage_counts, cell_types = stage_types, nUMI = stage_nUMI,
    min_UMI = 1, n_max_cells = Inf
  )
}
rm(chen_counts, chen_meta); gc()
log_mem("after building corrected references")

# ── Run RCTD for each of the 4 affected samples ─────────────────────────────
message("\n==> Running RCTD for the 4 affected samples ...")
weight_list <- vector("list", nrow(FIX_SAMPLES))
names(weight_list) <- FIX_SAMPLES$sample_id

for (i in seq_len(nrow(FIX_SAMPLES))) {
  row <- FIX_SAMPLES[i, ]
  message(sprintf("  [%d/%d] %s (stage %s)", i, nrow(FIX_SAMPLES), row$sample_id, row$stage))

  obj <- readRDS(file.path(DATA_DIR, row$rds_file))
  obj <- suppressWarnings(UpdateSeuratObject(obj))
  DefaultAssay(obj) <- "Spatial"
  obj[["SCT"]] <- NULL

  for (img_name in Images(obj)) {
    sf <- obj@images[[img_name]]@scale.factors$lowres
    obj@images[[img_name]]@coordinates$imagecol <- obj@images[[img_name]]@coordinates$imagecol * sf
    obj@images[[img_name]]@coordinates$imagerow <- obj@images[[img_name]]@coordinates$imagerow * sf
  }

  img <- obj@images[[1]]
  spatial_counts <- GetAssayData(obj, assay = "Spatial", layer = "counts")
  coords_df <- data.frame(x = img@coordinates$imagecol, y = img@coordinates$imagerow,
                           row.names = rownames(img@coordinates))
  spatial_nUMI <- colSums(spatial_counts)

  query <- SpatialRNA(coords_df, spatial_counts, spatial_nUMI)
  myRCTD <- create.RCTD(query, stage_refs[[row$stage]], max_cores = 16, CELL_MIN_INSTANCE = MIN_CELLS_PER_TYPE)
  rm(query); gc()
  myRCTD <- run.RCTD(myRCTD, doublet_mode = "full")
  w <- as.matrix(normalize_weights(myRCTD@results$weights))

  # Cell barcodes in the original merged object are prefixed with sample_id
  # (add.cell.ids at merge time in 05) -- reproduce that prefix so
  # AddMetaData matches the right cells below.
  rownames(w) <- paste0(row$sample_id, "_", rownames(w))

  w_full <- matrix(0, nrow = nrow(w), ncol = length(ALL_CELLTYPES),
                    dimnames = list(rownames(w), ALL_CELLTYPES))
  w_full[, colnames(w)] <- w
  weight_list[[row$sample_id]] <- w_full

  rm(myRCTD, obj, spatial_counts); gc()
  log_mem(paste("after RCTD for", row$sample_id))
}

# ── Load the deconvolved object and patch in the recovered weights ─────────
message("\n==> Loading VisiumLow_deconvolved.rds (patching in place) ...")
object <- readRDS(OBJECT_RDS)
log_mem("after loading deconvolved object")

fix_weights <- do.call(rbind, weight_list)
colnames(fix_weights) <- paste0("rctd_", colnames(fix_weights))
stopifnot(all(rownames(fix_weights) %in% colnames(object)))

fix_meta <- as.data.frame(fix_weights)
fix_meta$rctd_dominant_celltype <- sub(
  "^rctd_", "",
  colnames(fix_weights)[max.col(fix_weights, ties.method = "first")]
)

message(sprintf("  Patching %d cells (W12 + W52 samples)", nrow(fix_meta)))
object <- AddMetaData(object, fix_meta)

message("\n  Dominant cell type counts, W12/W52 cells, post-fix:")
print(table(object$rctd_dominant_celltype[rownames(fix_meta)], useNA = "ifany"))

message("\n==> Saving patched object ...")
saveRDS(object, OBJECT_RDS, compress = FALSE)
message("  Saved (overwritten): ", OBJECT_RDS)

message("\n==> Done.")
log_mem("final")
