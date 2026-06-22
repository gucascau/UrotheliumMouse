################################################################################
# 000_KidneySamples_MetaCorrection.R
#
# 
# Input files (both are scVI outputs — log-normalised expression in counts layer):
#   RenalUrothelium_uro_cells_fullgene_scvi.rds  —  11,105 kidney uro cells
#   
#
# Steps:
#  1. Load both scVI outputs as Seurat objects, ensuring data = counts (log-norm) and consistent metadata columns
#  2. Load the meta data: MetaDataCorrectionFinalVersion.csv
#  Ourput:
#  Save RenalUrothelium_uro_cells_fullgene_scvi_metacorrection.rds
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


OUT_PATH <- file.path(BASE_DIR, "RenalUrothelium_uro_cells_fullgene_scvi_metacorrection.rds")

################################################################################
# STEP 1: Load and set data = counts (scVI log-norm)
################################################################################

path     = file.path(BASE_DIR, "RenalUrothelium_uro_cells_fullgene_scvi.rds")
so <- readRDS(path)

# I need to add the meta data for the sample_id with Urotherlium 
gse_csv <- "/home/gdbecknelllab/xxw004/gdjacksonlab/UUO/Datasets/Mouse/UUO_CellMeta2022/GSE190887_meta_cell_type_sample.csv"

# Read the csv file
gse_meta <- read.csv(gse_csv, row.names = 1)
rownames(gse_meta) <- paste0(gsub("\\.", "_", rownames(gse_meta)), "-17")
gse_meta %>% head()
# add the cells
matched <- intersect(colnames(so), rownames(gse_meta))
so@meta.data$celltype0421     <- NA_character_
so@meta.data$gse190887_sample <- NA_character_
      # add the gsm_id

so@meta.data[matched, "celltype0421"]     <- gse_meta[matched, "celltype0421"]
so@meta.data[matched, "gse190887_sample"] <- gse_meta[matched, "sample"]

# check the condition column with the Urotherlium orig.ident

so@meta.data %>% filter(grepl("Urotherlium", orig.ident)) %>% select(condition) %>%table()
# I need to change the 

so@meta.data %>% filter(grepl("D2|D1", orig.ident)) %>% select(sample_id, orig.ident, condition, gsm_id) %>% head()

so@meta.data %>% filter(grepl("UUO", orig.ident)) %>% select(sample_id, orig.ident, condition,gsm_id) %>% tail()

so@meta.data %>% filter(orig.ident == "UUO") %>% select(sample_id, orig.ident, condition,gsm_id) %>% nrow()


# we change the orig.ident with the urothelum name and combine with condition
so@meta.data<- so@meta.data %>% mutate(orig.ident = case_when(
  # UUO_Day*/Healthy Urotherlium cells are already removed above; only IRI
  # Urotherlium cells remain — match all by name.
  # Note: the actual value is "Urotherlium" (with an extra 'r').
    orig.ident == "Urotherlium" ~ paste0(orig.ident, "_", gse190887_sample),
    orig.ident == "D1" ~ paste0(orig.ident, "_", condition),
    orig.ident == "D2" ~ paste0(orig.ident, "_", condition),
    orig.ident == "UUO" ~ "UUOKidney_14D",
   is.na(orig.ident) ~ sample_id,
  TRUE ~ orig.ident),
sample_id = case_when(
  grepl("Urotherlium", orig.ident) ~ orig.ident,
  grepl("D2|D1", orig.ident)  ~ orig.ident,
  grepl("UUO", orig.ident)  ~ orig.ident,
   is.na(sample_id) ~ orig.ident,
  TRUE ~ sample_id),
condition = case_when(
    grepl("Urotherlium", orig.ident) ~ gse190887_sample,
    grepl("D2|D1", orig.ident)  ~ condition,
    grepl("UUO", orig.ident)  ~ condition,
    TRUE ~ condition)
)



# check in the sample_id and orig.ident columns
so@meta.data$sample_id %>% table()
so@meta.data$orig.ident %>% table()

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
