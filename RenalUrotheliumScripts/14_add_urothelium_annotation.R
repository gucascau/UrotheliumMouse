#!/usr/bin/env Rscript
################################################################################
# 14_add_urothelium_annotation.R
#
# Add the urothelium subtype annotation (computed on the extracted urothelium
# object) back onto the full all-cells kidney object, matched by cell barcode.
#
# Input  : output/RenalUrothelium_allcells_scvi_annotations_metaupdated.rds
#          /vast0/.../FinalUrotheliumCells/RenalUrothelium_uro_cells_fullgene_scvi_metacorrection.rds
# Output : output/RenalUrothelium_allcells_scvi_annotations_meta_uro_updated.rds
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
OUT_DIR  <- file.path(BASE_DIR, "RenalUrotheliumScripts", "output")

IN_PATH   <- file.path(OUT_DIR, "RenalUrothelium_allcells_scvi_annotations_metaupdated.rds")
OUT_PATH  <- file.path(OUT_DIR, "RenalUrothelium_allcells_scvi_annotations_meta_uro_updated.rds")
URO_PATH  <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/FinalUrotheliumCells/RenalUrothelium_uro_cells_fullgene_scvi_metacorrection.rds"

# Urothelium-subtype columns computed only on the extracted urothelium object
# (not present anywhere in the all-cells object yet).
URO_COLS <- c("strict_urothelium", "uro_umbrella_arm", "uro_basal_arm",
              "uro_intermediate_arm", "uro_low_kidney_epi_score", "KidneyEpiScore1")

################################################################################
# STEP 1: Load both objects
################################################################################

cat("Reading", IN_PATH, "\n")
so <- readRDS(IN_PATH)
cat("  all-cells:", ncol(so), "cells\n")

cat("Reading", URO_PATH, "\n")
uro <- readRDS(URO_PATH)
cat("  urothelium:", ncol(uro), "cells\n")

missing_cols <- setdiff(URO_COLS, colnames(uro@meta.data))
if (length(missing_cols) > 0) {
  stop("Missing expected urothelium annotation column(s): ",
       paste(missing_cols, collapse = ", "))
}

################################################################################
# STEP 2: Match urothelium cells into the all-cells object by barcode
################################################################################

matched <- intersect(colnames(so), colnames(uro))
cat(sprintf("  %d / %d urothelium cells matched to all-cells barcodes\n",
            length(matched), ncol(uro)))
if (length(matched) < ncol(uro)) {
  cat("  WARNING: not all urothelium cells were found in the all-cells object.\n")
}

# Booleans default FALSE (not part of the strict urothelium call); scores stay NA.
so@meta.data$strict_urothelium        <- FALSE
so@meta.data$uro_umbrella_arm         <- FALSE
so@meta.data$uro_basal_arm            <- FALSE
so@meta.data$uro_intermediate_arm     <- FALSE
so@meta.data$uro_low_kidney_epi_score <- NA_real_
so@meta.data$KidneyEpiScore1          <- NA_real_

so@meta.data[matched, URO_COLS] <- uro@meta.data[matched, URO_COLS]

cat("\nstrict_urothelium table:\n")
print(table(so@meta.data$strict_urothelium, useNA = "ifany"))
cat("\nuro_umbrella_arm / uro_basal_arm / uro_intermediate_arm among strict_urothelium cells:\n")
print(so@meta.data %>%
  filter(strict_urothelium) %>%
  summarise(umbrella = sum(uro_umbrella_arm), basal = sum(uro_basal_arm),
            intermediate = sum(uro_intermediate_arm)))

################################################################################
# Save
################################################################################

cat("\nSaving to", OUT_PATH, "\n")
saveRDS(so, OUT_PATH)
cat("Done.\n")
