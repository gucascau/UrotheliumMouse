#!/usr/bin/env Rscript
################################################################################
# 13_allcells_metadata_update.R
#
# Correct sample_id/orig.ident/condition for the raw GSE190887 (IRI) samples
# ("Urotherlium", "D1", "D2", "UUO") the same way 000_KidneySamples_MetaCorrection.R
# does for the urothelium-only object, then left_join the harmonized Final*
# columns from MetaDataCorrectedFinalVersion.csv onto the full all-cells object.
#
# Input  : output/RenalUrothelium_allcells_scvi_annotations.rds
# Output : output/RenalUrothelium_allcells_scvi_annotations_metaupdated.rds
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

options(future.globals.maxSize = 8 * 1024^3)

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
OUT_DIR  <- file.path(BASE_DIR, "RenalUrotheliumScripts", "output")

IN_PATH  <- file.path(OUT_DIR, "RenalUrothelium_allcells_scvi_annotations.rds")
OUT_PATH <- file.path(OUT_DIR, "RenalUrothelium_allcells_scvi_annotations_metaupdated.rds")

GSE190887_CSV <- "/home/gdbecknelllab/xxw004/gdjacksonlab/UUO/Datasets/Mouse/UUO_CellMeta2022/GSE190887_meta_cell_type_sample.csv"
FINAL_META_CSV <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/FinalUrotheliumCells/UrotheliumScripts/output/AllUrothelium_metadata_correlations/MetaDataCorrectedFinalVersion_Update_07052026.csv"

################################################################################
# STEP 1: Load the all-cells object
################################################################################

cat("Reading", IN_PATH, "\n")
so <- readRDS(IN_PATH)
cat("  Cells:", ncol(so), "\n")


################################################################################
# STEP 2: Resolve the ambiguous raw sample_id groups (Urotherlium/D1/D2/UUO)
# using the GSE190887 per-cell sample lookup — same logic used for the
# urothelium-only object in 000_KidneySamples_MetaCorrection.R. The barcode
# suffix ("-17") matches because these cells were extracted from this same
# all-cells object, preserving the original batch index.
################################################################################

gse_meta <- read.csv(GSE190887_CSV, row.names = 1)
rownames(gse_meta) <- paste0(gsub("\\.", "_", rownames(gse_meta)), "-17")

matched <- intersect(colnames(so), rownames(gse_meta))
cat(sprintf("  GSE190887 lookup matched %d / %d Urotherlium-labeled cells\n",
            length(matched), sum(so@meta.data$orig.ident == "Urotherlium", na.rm = TRUE)))

so@meta.data$celltype0421     <- NA_character_
so@meta.data$gse190887_sample <- NA_character_
so@meta.data[matched, "celltype0421"]     <- gse_meta[matched, "celltype0421"]
so@meta.data[matched, "gse190887_sample"] <- gse_meta[matched, "sample"]

so@meta.data <- so@meta.data %>% mutate(orig.ident = case_when(
    # Note: the actual value is "Urotherlium" (with an extra 'r').
    orig.ident == "Urotherlium" ~ paste0(orig.ident, "_", gse190887_sample),
    orig.ident == "D1" ~ paste0(orig.ident, "_", condition),
    orig.ident == "D2" ~ paste0(orig.ident, "_", condition),
    orig.ident == "UUO" ~ "UUOKidney_14D",
    is.na(orig.ident) ~ sample_id,
    TRUE ~ orig.ident
  ),
  sample_id = case_when(
    grepl("Urotherlium", orig.ident) ~ orig.ident,
    grepl("D2|D1", orig.ident)  ~ orig.ident,
    grepl("UUO", orig.ident)  ~ orig.ident,
    is.na(sample_id) ~ orig.ident,
    TRUE ~ sample_id
  ),
  condition = case_when(
    grepl("Urotherlium", orig.ident) ~ gse190887_sample,
    grepl("D2|D1", orig.ident)  ~ condition,
    grepl("UUO", orig.ident)  ~ condition,
    TRUE ~ condition
  )
)

cat("  sample_id table after remapping:\n")
print(so@meta.data$sample_id %>% table())

# Check the metadata that we used for orig.ident and sample_id
cat("  sample_id table before remapping:\n")
print(so@meta.data$sample_id %>% table())
cat("  orig.ident table before remapping:\n")
print(so@meta.data$orig.ident %>% table())
so@meta.data$orig.ident %>% unique() %>% sort() %>% print()

# Before we applied the left_join of the Final* metadata, we need to remove the duplicated barcodes in the meta data, and keep the one with the correct sample_id, cells with sample_id "MouseKidneyATLAS_MKA_updated" are the duplicated one, we need to remove them from object, and keep the one with the correct sample_id
# there are some discrepancies in the sample_id and our meta data
so@meta.data %>% filter(grepl("MouseKidneyATLAS_MKA_updated", sample_id)) %>% head()

# check whether there are double ACCCACTCATTTGCTT
so@meta.data %>% filter(grepl("ACCCACTCATTTGCTT", rownames(so@meta.data))) %>% head()
so@meta.data %>% filter(grepl("TCGCGTTCAGCATGAG", rownames(so@meta.data))) %>% head()
# we need to remove the duplicated barcodes in the meta data, and keep the one with the correct sample_id, cells with sample_id "MouseKidneyATLAS_MKA_updated" are the duplicated one, we need to remove them from object, and keep the one with the correct sample_id
so <- so[, !grepl("MouseKidneyATLAS_MKA_updated", so@meta.data$sample_id)]


################################################################################
# STEP 3: Left-join the harmonized Final* metadata columns by sample_id
################################################################################

UpdateMetaData <- read.csv(FINAL_META_CSV)

new_cols <- setdiff(colnames(UpdateMetaData), colnames(so@meta.data))
new_meta <- so@meta.data %>%
  left_join(select(UpdateMetaData, all_of(c("sample_id", new_cols))), by = "sample_id")
rownames(new_meta) <- rownames(so@meta.data)
so@meta.data <- new_meta

################################################################################
# STEP 4: Diagnostics — sample_id groups with no match in the Final metadata
# (left as NA; not treated as an error — see conversation notes).
################################################################################

unmatched <- so@meta.data %>%
  filter(is.na(FinalConditionL1)) %>%
  count(sample_id, sort = TRUE)

cat("\nSample_id groups with NO match in", basename(FINAL_META_CSV), "(Final* columns left NA):\n")
print(unmatched)

cat("\nsample_id values present in object but absent from the Final metadata CSV:\n")
print(setdiff(unique(so@meta.data$sample_id), unique(UpdateMetaData$sample_id)))

################################################################################
# STEP 5: Check the Final* columns for NA values (should be none, except for the unmatched sample_id groups above)
# Check the final annotation columns for NA values (should be none, except for the unmatched sample_id groups above)
################################################################################
# Check the Final* columns for NA values (should be none, except for the unmatched sample_id groups above)
library(tidyr)
una_final_cols <- so@meta.data %>%
  select(starts_with("Final")) %>%
  summarise_all(~ sum(is.na(.))) %>%
  pivot_longer(everything(), names_to = "column", values_to = "num_na") %>%
  filter(num_na > 0)

so@meta.data %>% filter(is.na(Finalgsm_id)) %>% count(sample_id, sort = TRUE) %>% print()
# Check the annotation columns for NA values (should be none, except for the unmatched sample_id groups above), the column name contain annotation, so we can use the grepl to filter the column name
colnames(so@meta.data)
so@meta.data %>% head()
so@meta.data %>% select(matches("celltype|MKA|Lake|cell_type|suspension_type|scanvi")) %>% head()

# print out the MKA and Lake annotation numbers
so@meta.data %>% count(scanvi_MKA_label, sort = TRUE) %>% print()
so@meta.data %>% count(scanvi_Lake_label, sort = TRUE) %>% print()

# we finally used MKA label as the final annotation, so we can check the MKA label for NA values
so@meta.data %>% filter(is.na(scanvi_MKA_label)) %>% count(sample_id, sort = TRUE) %>% print()


unannotated_cols <- so@meta.data %>%
  select(where(~ grepl("annotation", .))) %>%
  summarise_all(~ sum(is.na(.))) %>%
  pivot_longer(everything(), names_to = "column", values_to = "num_na") %>%
  filter(num_na > 0)

################################################################################
# Save
################################################################################

cat("\nSaving to", OUT_PATH, "\n")
saveRDS(so, OUT_PATH)
cat("Done.\n")
