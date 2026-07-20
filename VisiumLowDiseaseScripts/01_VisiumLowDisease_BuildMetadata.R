################################################################################
# 01_VisiumLowDisease_BuildMetadata.R
#
# Build authoritative per-sample metadata for the 10 standard-resolution
# ("low-res", 55um-spot) adult mouse kidney Visium samples living alongside
# the developmental atlas in UsedSpatialData/VisiumLow/, but NOT part of it --
# these are disease-model and healthy-adult-reference datasets from three
# different sources, kept as a deliberately separate integration from
# VisiumLowDevelopmentScripts (embryonic->aged time course, no injury):
#
#   - GSE269063 (2 samples): UUO surgery study, Sham vs Vehicle. Raw GEO
#     supplementary files (no spaceranger outs/ folder) -- whole-transcriptome
#     chemistry, mm10-2020-A reference (32,285 genes).
#   - GSE269884 (6 samples): bilateral IRI (Bi-IRI) injury time course --
#     sham, 4h, 12h, 2d, 14d, 6wk post-surgery, all male. Full spaceranger
#     outs/ directories. 10X Visium CytAssist, Visium Mouse Transcriptome
#     Probe Set v1.0 mm10-2020-A (19,465 genes; PROBE-BASED, not
#     whole-transcriptome).
#   - 10x Genomics public Datasets (2 samples, not GEO): "V1_Mouse_Kidney"
#     (Visium_Coranal_Mouse_Kidney/, fresh-frozen coronal section, healthy
#     adult, whole-transcriptome, 32,285 genes) and "Visium_FFPE_Mouse_Kidney"
#     (FFPE, healthy adult, same probe panel as GSE269884, 19,465 genes).
#     No GSM/soft file to parse -- provenance is the 10x product page file
#     naming only, recorded here rather than guessed.
#
# Confirmed empirically (see conversation) that the probe panel (19,465
# genes) is an exact subset of the whole-transcriptome gene set (32,285
# genes) in both probe/whole-transcriptome pairs -- same mm10-2020-A build
# throughout. 02_VisiumLowDisease_Integrate_Harmony.R intersects every sample
# down to the 19,465-gene probe panel before merging, since Harmony cannot
# meaningfully batch-correct across genes present in some samples and absent
# (structurally, not just undetected) in others.
#
# load_type dispatch (used by 02_VisiumLowDisease_Integrate_Harmony.R):
#   - "raw_geo_flat"     : flat, gzipped GEO supplementary files, no outs/
#                          folder -- matrix/barcodes/features + a separate
#                          spatial/ set of files to stage. (GSE269063)
#   - "spaceranger_outs" : a complete spaceranger outs/ directory already on
#                          disk -- filtered_feature_bc_matrix.h5 + spatial/.
#                          (GSE269884)
#   - "h5_plus_spatial_tar": a filtered_feature_bc_matrix.h5 already extracted
#                          on disk, but spatial/ still packed in a
#                          spatial.tar.gz needing extraction. (10x Datasets)
#
# Output: VisiumLowDiseaseScripts/output/VisiumLowDisease_sample_metadata.csv
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

DATA_DIR      <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSpatialData/VisiumLow"
GSE269063_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/RenalDiseases/Mouse/UUO/GSE269063"
SCRIPT_DIR    <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/VisiumLowDiseaseScripts"
OUT_DIR       <- file.path(SCRIPT_DIR, "output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Parse GEO soft files for the two GEO series (traceable, not hand-typed) ──
parse_soft_chars <- function(soft_gz, gsm_ids) {
  lines <- readLines(gzfile(soft_gz))
  start_idx <- grep("^\\^SAMPLE = ", lines)
  end_idx   <- c(start_idx[-1] - 1, length(lines))
  acc_all   <- str_match(lines[start_idx], "\\^SAMPLE = (GSM[0-9]+)")[, 2]
  keep      <- which(acc_all %in% gsm_ids)

  do.call(rbind, lapply(keep, function(i) {
    block <- lines[start_idx[i]:end_idx[i]]
    title <- str_match(block, "^!Sample_title = (.*)$")
    title <- title[!is.na(title[, 1]), 2][1]
    chars <- block[grepl("^!Sample_characteristics_ch1 = ", block)]
    chars <- str_remove(chars, "^!Sample_characteristics_ch1 = ")
    get_char <- function(key) {
      hit <- chars[str_starts(chars, fixed(paste0(key, ":")))]
      if (length(hit) == 0) return(NA_character_)
      str_trim(str_remove(hit[1], paste0("^", key, ":")))
    }
    data.frame(
      GSM          = acc_all[i],
      title        = title,
      tissue       = get_char("tissue"),
      condition    = get_char("condition"),
      treatment    = get_char("treatment"),
      sex_raw      = get_char("Sex"),
      genotype     = get_char("genotype"),
      stringsAsFactors = FALSE
    )
  }))
}

message("==> Parsing GSE269063_family.soft.gz (UUO, 2 samples) ...")
uuo_meta <- parse_soft_chars(
  file.path(GSE269063_DIR, "GSE269063_family.soft.gz"),
  c("GSM8305769", "GSM8305770")
)
stopifnot(nrow(uuo_meta) == 2)

message("==> Parsing GSE269884_family.soft.gz (Bi-IRI time course, 6 samples) ...")
iri_meta <- parse_soft_chars(
  file.path(DATA_DIR, "..", "GSE269884_family.soft.gz"),
  sprintf("GSM%d", 8323120:8323125)
)
stopifnot(nrow(iri_meta) == 6)

# ── Assemble the 10-sample table ─────────────────────────────────────────────
# GSE269063: raw flat GEO supplementary files, whole-transcriptome
uuo_meta <- uuo_meta %>%
  mutate(
    sample_id     = ifelse(GSM == "GSM8305769", "GSE269063_SH045", "GSE269063_Vehicle"),
    study         = "GSE269063",
    disease_model = "UUO",
    timepoint     = treatment,
    timepoint_hours = NA_real_,
    sex           = NA_character_,  # not reported in this series' characteristics
    chemistry     = "WholeTranscriptome",
    tissue_prep   = "FreshFrozen",
    n_genes_expected = 32285L,
    load_type     = "raw_geo_flat",
    matrix_dir    = GSE269063_DIR,
    matrix_prefix = ifelse(GSM == "GSM8305769", "GSM8305769_003_SH045", "GSM8305770_002_Veh"),
    spatial_dir   = NA_character_,
    outs_dir      = NA_character_,
    h5_path       = NA_character_,
    spatial_tar   = NA_character_
  )

# GSE269884: full spaceranger outs/ dirs already on disk, CytAssist probe-based
iri_outs_dirname <- c(
  GSM8323120 = "GSM8323120_visium_sham_male",
  GSM8323121 = "GSM8323121_visium_hour4_male",
  GSM8323122 = "GSM8323122_visium_hour12_male",
  GSM8323123 = "GSM8323123_visium_day2_male",
  GSM8323124 = "GSM8323124_visium_day14_male",
  GSM8323125 = "GSM8323125_visium_week6_male"
)
iri_timepoint_hours <- c(
  GSM8323120 = 0,     # sham, day14-matched contemporaneous control
  GSM8323121 = 4,
  GSM8323122 = 12,
  GSM8323123 = 48,
  GSM8323124 = 336,
  GSM8323125 = 1008
)
iri_meta <- iri_meta %>%
  mutate(
    sample_id     = paste0("IRI_", str_remove(title, "_male$")),
    study         = "GSE269884",
    disease_model = "Bilateral_IRI",
    timepoint     = str_remove(title, "_male$"),
    timepoint_hours = iri_timepoint_hours[GSM],
    sex           = "male",
    chemistry     = "ProbeBasedCytAssist",
    tissue_prep   = "FreshFrozen",
    n_genes_expected = 19465L,
    load_type     = "spaceranger_outs",
    matrix_dir    = NA_character_,
    matrix_prefix = NA_character_,
    spatial_dir   = NA_character_,
    outs_dir      = file.path(DATA_DIR, iri_outs_dirname[GSM], "outs"),
    h5_path       = NA_character_,
    spatial_tar   = NA_character_
  )

# 10x Genomics public reference datasets (no GEO record -- provenance is the
# product-page file naming, recorded directly rather than inferred).
refs_meta <- data.frame(
  GSM              = NA_character_,
  title             = c("V1_Mouse_Kidney (10x Genomics public dataset)",
                         "Visium_FFPE_Mouse_Kidney (10x Genomics public dataset)"),
  tissue            = "kidney",
  condition         = "Healthy_Reference",
  treatment          = NA_character_,
  sex_raw           = NA_character_,
  genotype          = NA_character_,
  sample_id         = c("Coranal_HealthyRef", "FFPE_HealthyRef"),
  study             = "10x_Genomics_Datasets",
  disease_model     = "None",
  timepoint         = NA_character_,
  timepoint_hours   = NA_real_,
  sex               = NA_character_,
  chemistry         = c("WholeTranscriptome", "ProbeBasedFFPE"),
  tissue_prep       = c("FreshFrozen", "FFPE"),
  n_genes_expected  = c(32285L, 19465L),
  load_type         = "h5_plus_spatial_tar",
  matrix_dir        = NA_character_,
  matrix_prefix     = NA_character_,
  spatial_dir       = NA_character_,
  outs_dir          = NA_character_,
  h5_path           = c(
    file.path(DATA_DIR, "Visium_Coranal_Mouse_Kidney", "V1_Mouse_Kidney_filtered_feature_bc_matrix.h5"),
    file.path(DATA_DIR, "Visium_FFPE_Mouse_Kidney", "Visium_FFPE_Mouse_Kidney_filtered_feature_bc_matrix.h5")
  ),
  spatial_tar       = c(
    file.path(DATA_DIR, "Visium_Coranal_Mouse_Kidney", "V1_Mouse_Kidney_spatial.tar.gz"),
    file.path(DATA_DIR, "Visium_FFPE_Mouse_Kidney", "Visium_FFPE_Mouse_Kidney_spatial.tar.gz")
  ),
  stringsAsFactors = FALSE
)

sample_meta <- bind_rows(uuo_meta, iri_meta, refs_meta) %>%
  select(sample_id, GSM, study, disease_model, condition, timepoint, timepoint_hours,
         sex, chemistry, tissue_prep, n_genes_expected, load_type,
         matrix_dir, matrix_prefix, spatial_dir, outs_dir, h5_path, spatial_tar,
         title, treatment, genotype)

stopifnot(nrow(sample_meta) == 10)
stopifnot(!anyDuplicated(sample_meta$sample_id))

message("\n  Final metadata table:")
print(sample_meta %>% select(sample_id, study, disease_model, timepoint, chemistry, load_type))

message("\n  Samples per study:")
print(table(sample_meta$study))
message("\n  Samples per chemistry:")
print(table(sample_meta$chemistry))

out_csv <- file.path(OUT_DIR, "VisiumLowDisease_sample_metadata.csv")
write.csv(sample_meta, out_csv, row.names = FALSE)
message("\n  Saved: ", out_csv)
