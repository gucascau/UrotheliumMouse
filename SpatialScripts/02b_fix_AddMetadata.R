################################################################################
# 02b_fix_AddMetadata.R
#
# Recovery script for 02b_VisiumHD_Deconvolution_BothKidneys.R.
#
# The RCTD loops completed successfully and saved a checkpoint:
#   output/all_rctd_meta_bothkidneys.rds
#
# The job then failed at AddMetaData() because RCTD "full" mode silently drops
# low-UMI bins, so nrow(combined_meta) < ncol(object) and AddMetaData rejects
# a data frame that does not cover every cell.
#
# Fix: index into combined_meta by cell barcode (rownames of @meta.data).
# Bins absent from RCTD results receive NA — correct behaviour.
#
# Outputs (same as the original script):
#   output/VisiumHD_kidney_both_deconvolved.rds
#   output/VisiumHD_RCTD_dominant_celltype.pdf
################################################################################

library(Seurat)
library(ggplot2)
library(dplyr)

options(future.globals.maxSize = 16 * 1024^3)

log_mem <- function(label) {
  g  <- gc(verbose = FALSE)
  mb <- sum(g[, 2])
  message(sprintf("  [mem] %s: %.1f GB in use", label, mb / 1024))
}

OUT_DIR    <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/SpatialScripts/output"
OBJECT_RDS <- file.path(OUT_DIR, "VisiumHD_harmony_KidneyOnlyintegrated.rds")
CKPT_RDS   <- file.path(OUT_DIR, "all_rctd_meta.rds")
OUT_DECONV <- file.path(OUT_DIR, "VisiumHD_kidney_both_deconvolved.rds")
OUT_PDF    <- file.path(OUT_DIR, "VisiumHD_RCTD_dominant_celltype.pdf")

# we have two side of the kidneys that have been annotated 
# Left kidney
LeftKidneyAnnObj <- file.path(OUT_DIR, "VisiumHD_kidney3p_urothelium_deconvolved.rds")

# Right Kidney
RightKidneyAnnObj <- file.path(OUT_DIR, "VisiumHD_kidney5p_deconvolved.rds")


if (!file.exists(CKPT_RDS))
  stop("Checkpoint not found: ", CKPT_RDS,
       "\nRe-run 02b_VisiumHD_Deconvolution_BothKidneys.R to regenerate RCTD results.")

# ── Load ──────────────────────────────────────────────────────────────────────
message("==> Loading VisiumHD integrated object ...")
object <- readRDS(OBJECT_RDS)
DefaultAssay(object) <- "Spatial.008um"
log_mem("after loading object")
# check the meta data for the object
object %>% head()


message("==> Loading VisiumHD left and right kidney objects ...")
LeftKidneyAnnObject <- readRDS(LeftKidneyAnnObj)
RightKidneyAnnObject <- readRDS(RightKidneyAnnObj)

LeftKidneyAnnObject %>% head()
RightKidneyAnnObject %>% head()

LeftKidneyAnnObject@meta.data %>% pull(rctd_dominant_celltype) %>% table()

# check how many of the cell with NA for rctd_dominant_celltype
LeftKidneyAnnObject@meta.data %>% filter(is.na(rctd_dominant_celltype)) %>% nrow()
LeftKidneyAnnObject@meta.data %>% filter(!is.na(rctd_dominant_celltype)) %>% nrow()
# check how many fo the cell with NA in 
RightKidneyAnnObject@meta.data %>% filter(is.na(rctd_dominant_celltype)) %>% nrow()
RightKidneyAnnObject@meta.data %>% filter(!is.na(rctd_dominant_celltype)) %>% nrow()

# we save the cell id and annotated cell types 
LeftKidneyCellAnnotation <- LeftKidneyAnnObject@meta.data %>% select()

log_mem("after loading object")

message("==> Loading RCTD checkpoint ...")
all_rctd_meta <- readRDS(CKPT_RDS)
message(sprintf("  Images in checkpoint: %s", paste(names(all_rctd_meta), collapse = ", ")))
for (img in names(all_rctd_meta)) {
  message(sprintf("  %-30s: %d bins with RCTD weights", img, nrow(all_rctd_meta[[img]])))
}

# ── Combine RCTD results ───────────────────────────────────────────────────────
message("\n==> Combining RCTD results ...")
combined_meta <- do.call(rbind, unname(all_rctd_meta))
rm(all_rctd_meta); gc()

message(sprintf("  combined_meta: %d rows,  object cells: %d",
                nrow(combined_meta), ncol(object)))

n_missing <- sum(!colnames(object) %in% rownames(combined_meta))
if (n_missing > 0)
  message(sprintf("  %d bins have no RCTD weights (filtered by RCTD) — will be NA",
                  n_missing))

# ── Add metadata safely by barcode matching ───────────────────────────────────
# Direct index on rownames(object@meta.data) → NA for bins absent from RCTD.
message("\n==> Adding RCTD weights to Seurat object ...")
for (col in colnames(combined_meta)) {
  object@meta.data[[col]] <- combined_meta[rownames(object@meta.data), col]
}
rm(combined_meta); gc()
log_mem("after adding metadata")

message("Dominant cell type distribution:")
print(table(object$rctd_dominant_celltype, useNA = "ifany"))

# ── Save ──────────────────────────────────────────────────────────────────────
message("\n==> Saving deconvolved object ...")
saveRDS(object, OUT_DECONV)
message("  Saved: ", OUT_DECONV)

# ── Spatial plots — one panel per kidney sample ───────────────────────────────
message("==> Plotting ...")
all_images <- Images(object)
pdf(OUT_PDF, width = 10, height = 7)

for (img in all_images) {
  tryCatch({
    p <- SpatialDimPlot(object,
      group.by       = "rctd_dominant_celltype",
      images         = img,
      crop           = TRUE,
      pt.size.factor = 2,
      image.alpha    = 1
    ) + ggtitle(paste("RCTD cell types:", img)) +
      theme(legend.position = "right")
    print(p)
  }, error = function(e) {
    message("  SpatialDimPlot failed for ", img, ": ", e$message)
  })
}

dev.off()
message("  Saved: ", OUT_PDF)

message("\n==> Done.")
log_mem("final")
