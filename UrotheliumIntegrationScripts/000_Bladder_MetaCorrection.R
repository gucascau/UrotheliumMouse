################################################################################
# 000_Bladderamples_MetaCorrection.R
#
# 
# Input files (both are scVI outputs — log-normalised expression in counts layer):
#   BladderUrothelium_uro_cells_scvi.rds  
#   
#
# Steps:
#  1. Load both scVI outputs as Seurat objects, ensuring data = counts (log-norm) and consistent metadata columns
#  2. Load the meta data: MetaDataCorrectionFinalVersion.csv
#  Ourput:
#  Save BladderUrothelium_uro_cells_scvi_metacorrection.rds
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


OUT_PATH <- file.path(BASE_DIR, "BladderUrothelium_uro_cells_scvi_metacorrection.rds")

################################################################################
# STEP 1: Load and set data = counts (scVI log-norm)
################################################################################

path     = file.path(BASE_DIR, "BladderUrothelium_uro_cells_scvi.rds")
so <- readRDS(path)
so@meta.data %>% head()

# check the sample_id
so@meta.data$sample_id %>% table()
# check the BladderNormal1, BladderNormal1
so@meta.data %>% filter(sample_id %in% c("BladderNormal1", "BladderNormal2")) %>% head()

# we read the meta data file and concatenate by the sample_id
UpdateMetaData <- read.csv("/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/FinalUrotheliumCells/UrotheliumScripts/output/AllUrothelium_metadata_correlations/MetaDataCorrectedFinalVersion.csv")

# Only bring in columns that are new (avoid duplicate column conflicts).
# left_join drops rownames, so restore them from the original meta.data.
UpdateMetaData$sample_id %>% table()

new_cols <- setdiff(colnames(UpdateMetaData), colnames(so@meta.data))
new_meta <- so@meta.data %>%
  left_join(select(UpdateMetaData, all_of(c("sample_id", new_cols))), by = "sample_id")
rownames(new_meta) <- rownames(so@meta.data)
so@meta.data <- new_meta
setdiff(unique(so@meta.data$sample_id), unique(UpdateMetaData$sample_id))

new_meta %>% filter(is.na(FinalConditionL1)) %>% tail()
so@meta.data %>% select(sample_id, orig.ident, condition, tissue, technology, paper, FinalConditionL2, FinalConditionL1, Finalgsm_id, Finalpaper ) %>% tail()

# Save the updated Seurat object with corrected metadata
saveRDS(so, OUT_PATH)
