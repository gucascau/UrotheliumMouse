################################################################################
# 04_merge_GEO_metadata.R
#
# Left-join AllUrothelium_orig_ident_GEO_metadata.csv (GEO accessions, one row
# per orig.ident) with MetaDataNeedCorretion.csv (condition levels, categories,
# sample_id) on orig.ident.
#
# Note: "Urotherlium" maps to 17 rows in MetaDataNeedCorretion (one per
# condition/timepoint in the GSE190887 KUDO dataset), so the output has
# 123 rows rather than 107.
#
# Output: output/AllUrothelium_orig_ident_GEO_metadata_merged.csv
################################################################################

library(dplyr)
library(readr)

OUT_DIR <- "output"

# ── Load ──────────────────────────────────────────────────────────────────────

geo  <- read_csv(file.path(OUT_DIR, "AllUrothelium_orig_ident_GEO_metadata.csv"),
                 show_col_types = FALSE)

corr <- read_csv(file.path(OUT_DIR, "MetaDataNeedCorretion.csv"),
                 show_col_types = FALSE)

# Drop any trailing empty columns introduced by Excel/LibreOffice
corr <- corr |> select(where(~ !all(is.na(.))))

message(sprintf("GEO metadata  : %d rows × %d cols", nrow(geo),  ncol(geo)))
message(sprintf("Correction CSV: %d rows × %d cols", nrow(corr), ncol(corr)))

# ── Check overlap ─────────────────────────────────────────────────────────────

only_in_geo  <- setdiff(geo$orig.ident,  corr$orig.ident)
only_in_corr <- setdiff(corr$orig.ident, geo$orig.ident)

if (length(only_in_geo)  > 0) message("In GEO only : ", paste(only_in_geo,  collapse = ", "))
if (length(only_in_corr) > 0) message("In Corr only: ", paste(only_in_corr, collapse = ", "))

# ── Left join ─────────────────────────────────────────────────────────────────
# Columns kept unchanged from AllUrothelium_orig_ident_GEO_metadata.csv:
#   source_GEO, gsm_id, paper, technology, tissue, condition, notes
#
# Overlapping columns from MetaDataNeedCorretion.csv are added as NEW columns
# with a "_corr" suffix (not replacing the originals):
#   source_GEO_corr, gsm_id_corr, paper_corr, technology_corr,
#   tissue_corr, condition_corr
#
# Columns added from MetaDataNeedCorretion.csv only:
#   sample_id, Categories, condition_level1/2/3, FinalConditionL1/2

# Rename overlapping columns in corr before joining
overlap_cols <- c("source_GEO", "gsm_id", "paper", "technology", "tissue", "condition")

corr_renamed <- corr |>
  rename_with(~ paste0(.x, "_corr"), .cols = all_of(overlap_cols))

merged <- geo |>
  left_join(corr_renamed, by = "orig.ident")

message(sprintf("Merged result : %d rows × %d cols", nrow(merged), ncol(merged)))

# Summarise any orig.ident that expanded (1-to-many)
expansion <- merged |>
  count(orig.ident, name = "n_rows") |>
  filter(n_rows > 1)

if (nrow(expansion) > 0) {
  message("\nExpanded orig.ident values (multiple rows in MetaDataNeedCorretion):")
  print(expansion, n = Inf)
}

# ── Save ──────────────────────────────────────────────────────────────────────

out_path <- file.path(OUT_DIR, "AllUrothelium_orig_ident_GEO_metadata_merged.csv")
write_csv(merged, out_path)
message(sprintf("\nSaved → %s", out_path))
