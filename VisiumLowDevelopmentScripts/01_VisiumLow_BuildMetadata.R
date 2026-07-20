################################################################################
# 01_VisiumLow_BuildMetadata.R
#
# Build authoritative per-sample metadata for the 24 standard-resolution
# ("low-res", 55um-spot) Visium developmental-kidney samples (GSE252772,
# GSM8704836-GSM8704859) -- the spatial companion to the Chen2025 NatGenet
# snRNA/snATAC developmental atlas used throughout UrotheliumDevelopmentScripts
# (same study, same 6 stages: E16.5, P0, W3, W12, W52, W92).
#
# Parsed directly from the GEO family.soft.gz (not hand-typed / not guessed
# from filenames) so the metadata is reproducible and traceable to source.
# Filenames already carry the GSM accession as a prefix
# ("GSM8704836_20201201-..._obj.rds"), which is matched exactly against the
# soft file's ^SAMPLE blocks -- no fuzzy string matching.
#
# Two known complications, confirmed from the GEO Sample_description field
# (see README note below and memory):
#   - GSM8704846, 849, 850, 851, 852 (all nominally "P0"): each Visium slide
#     actually contains a female kidney section AND a male kidney section
#     side by side ("Female.Left"/"Male.Right" in the description). A single
#     per-sample Sex value is wrong for these -- flagged as sex = "Mixed" and
#     mixed_sex_slide = TRUE rather than guessed. Any sex-stratified analysis
#     should exclude or specially handle these 5 samples.
#   - GSM8704855, 856, 857 (all "E16.5"): each slide pools 4 separate embryo
#     sections (one animal, same sex throughout that slide) rather than one
#     specimen -- captured in n_sections_pooled, not treated as 4 independent
#     biological replicates.
#
# Input:  UsedSpatialData/VisiumLow/DevelopedKidney/GSE252772_family.soft.gz
#         UsedSpatialData/VisiumLow/DevelopedKidney/*_obj.rds (to confirm
#         every RDS file has a matching metadata row, and vice versa)
# Output: VisiumLowDevelopmentScripts/output/VisiumLow_sample_metadata.csv
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

DATA_DIR   <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSpatialData/VisiumLow/DevelopedKidney"
SCRIPT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/VisiumLowDevelopmentScripts"
OUT_DIR    <- file.path(SCRIPT_DIR, "output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

SOFT_GZ <- file.path(DATA_DIR, "GSE252772_family.soft.gz")

STAGE_ORDER <- c("E16.5", "P0", "W3", "W12", "W52", "W92")
AGE_MAP <- c(E16 = "E16.5", P0 = "P0", W3 = "W3", W12 = "W12", W52 = "W52", W92 = "W92")

# Visium slides with two different-sex specimens pooled on one capture area
# (per GEO Sample_description "...Left); ...(...Right)" -- confirmed by hand
# above, hardcoded here since it's a fact about the source data, not derived).
MIXED_SEX_GSM <- c("GSM8704846", "GSM8704849", "GSM8704850", "GSM8704851", "GSM8704852")

# ── Parse the GEO soft file ──────────────────────────────────────────────────
message("==> Parsing ", basename(SOFT_GZ), " ...")
lines <- readLines(gzfile(SOFT_GZ))

sample_start_idx <- grep("^\\^SAMPLE = ", lines)
sample_end_idx    <- c(sample_start_idx[-1] - 1, length(lines))

parse_block <- function(block) {
  acc <- str_match(block[1], "\\^SAMPLE = (GSM[0-9]+)")[, 2]
  title <- str_match(block, "^!Sample_title = (.*)$")
  title <- title[!is.na(title[, 1]), 2][1]
  chars <- block[grepl("^!Sample_characteristics_ch1 = ", block)]
  chars <- str_remove(chars, "^!Sample_characteristics_ch1 = ")
  descs <- block[grepl("^!Sample_description = ", block)]
  descs <- str_remove(descs, "^!Sample_description = ")

  get_char <- function(key) {
    hit <- chars[str_starts(chars, fixed(paste0(key, ":")))]
    if (length(hit) == 0) return(NA_character_)
    str_trim(str_remove(hit[1], paste0("^", key, ":")))
  }

  library_name <- str_match(descs[1], "Library name: (.*)$")[, 2]
  specimen_desc <- if (length(descs) >= 2) descs[2] else NA_character_
  n_pooled <- if (!is.na(specimen_desc)) length(str_split(specimen_desc, ";\\s*")[[1]]) else 1L

  data.frame(
    GSM = acc,
    title = title,
    library_name = library_name,
    tissue = get_char("tissue"),
    age_raw = get_char("age"),
    sex_raw = get_char("Sex"),
    strain = get_char("strain"),
    genotype = get_char("genotype"),
    specimen_description = specimen_desc,
    n_sections_pooled = n_pooled,
    stringsAsFactors = FALSE
  )
}

all_samples <- do.call(rbind, lapply(seq_along(sample_start_idx), function(i) {
  parse_block(lines[sample_start_idx[i]:sample_end_idx[i]])
}))

# ── Restrict to the 24 Visium samples (GSM8704836-859) ──────────────────────
visium_meta <- all_samples %>% filter(GSM %in% sprintf("GSM%d", 8704836:8704859))
message(sprintf("  Parsed %d Visium sample records (expected 24)", nrow(visium_meta)))
stopifnot(nrow(visium_meta) == 24)

# ── Normalize Age to the same labels used throughout UrotheliumDevelopmentScripts ──
# GEO uses "W3"/"E16"/etc directly for the two age-format characteristics rows
# for these Visium samples (unlike the paired snRNA samples, which had a
# separate "3w"/"e16.5" row too) -- map through AGE_MAP for E16 -> E16.5,
# identity otherwise.
visium_meta$Age <- ifelse(visium_meta$age_raw %in% names(AGE_MAP),
                           AGE_MAP[visium_meta$age_raw], visium_meta$age_raw)
stopifnot(all(visium_meta$Age %in% STAGE_ORDER))
visium_meta$Age <- factor(visium_meta$Age, levels = STAGE_ORDER)

# ── Sex, with the mixed-slide correction ─────────────────────────────────────
visium_meta$mixed_sex_slide <- visium_meta$GSM %in% MIXED_SEX_GSM
visium_meta$Sex <- ifelse(visium_meta$mixed_sex_slide, "Mixed", visium_meta$sex_raw)

# ── Match to the actual RDS filenames on disk ────────────────────────────────
rds_files <- list.files(DATA_DIR, pattern = "_obj\\.rds$", full.names = FALSE)
rds_gsm   <- str_match(rds_files, "^(GSM[0-9]+)_")[, 2]
file_map  <- setNames(rds_files, rds_gsm)

missing_rds <- setdiff(visium_meta$GSM, names(file_map))
missing_meta <- setdiff(names(file_map), visium_meta$GSM)
if (length(missing_rds) > 0) stop("No RDS file found for GSM(s): ", paste(missing_rds, collapse = ", "))
if (length(missing_meta) > 0) stop("No metadata found for RDS file GSM(s): ", paste(missing_meta, collapse = ", "))

visium_meta$rds_file <- file_map[visium_meta$GSM]
visium_meta$sample_id <- str_remove(visium_meta$rds_file, "_obj\\.rds$")

visium_meta <- visium_meta %>%
  select(GSM, sample_id, rds_file, Age, Sex, mixed_sex_slide, n_sections_pooled,
         title, library_name, specimen_description, tissue, strain, genotype) %>%
  arrange(Age, GSM)

message("\n  Final metadata table:")
print(visium_meta %>% select(GSM, sample_id, Age, Sex, mixed_sex_slide, n_sections_pooled))

message("\n  Samples per stage:")
print(table(visium_meta$Age))
message("\n  Mixed-sex slides: ", sum(visium_meta$mixed_sex_slide), " (as expected: 5)")

out_csv <- file.path(OUT_DIR, "VisiumLow_sample_metadata.csv")
write.csv(visium_meta, out_csv, row.names = FALSE)
message("\n  Saved: ", out_csv)
