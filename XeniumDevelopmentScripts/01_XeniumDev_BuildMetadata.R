################################################################################
# 01_XeniumDev_BuildMetadata.R
#
# Build per-sample metadata for the 4 Xenium mouse-kidney samples (GSE286051)
# -- W12 and W92, one male + one female each. Same overall study as
# VisiumLowDevelopmentScripts/UrotheliumDevelopmentScripts (Chen2025 NatGenet,
# "sex-specific differences in gene regulation across the lifespan"), parsed
# from GSE286051_family.soft.gz -- not the 12-sample UUO time-course Xenium
# dataset already handled by SpatialScripts/02_Xenium_integrate_harmony.R.
#
# Files are flat GEO supplementary files (matrix.mtx.gz/barcodes.tsv.gz/
# features.tsv.gz/cells.csv.gz/transcripts.parquet.gz/morphology.ome.tif.gz),
# NOT a standard Xenium outs/ bundle -- no cell_boundaries.parquet or
# experiment.xenium, so LoadXenium() can't be used directly.
# cells.csv.gz only gives centroid coordinates (x_centroid/y_centroid), no
# polygon segmentation -- 02_XeniumDev_Integrate_Harmony.R builds a
# centroids-only FOV (CreateCentroids + CreateFOV(type = "centroids")).
#
# Input:  UsedSpatialData/Xenium/GSE286051/GSE286051_family.soft.gz
# Output: XeniumDevelopmentScripts/output/XeniumDev_sample_metadata.csv
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

DATA_DIR   <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSpatialData/Xenium/GSE286051"
SCRIPT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/XeniumDevelopmentScripts"
OUT_DIR    <- file.path(SCRIPT_DIR, "output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

SOFT_GZ <- file.path(DATA_DIR, "GSE286051_family.soft.gz")
GSM_IDS <- c("GSM8717006", "GSM8717007", "GSM8717008", "GSM8717009")

message("==> Parsing ", basename(SOFT_GZ), " ...")
lines <- readLines(gzfile(SOFT_GZ))
start_idx <- grep("^\\^SAMPLE = ", lines)
end_idx   <- c(start_idx[-1] - 1, length(lines))
acc_all   <- str_match(lines[start_idx], "\\^SAMPLE = (GSM[0-9]+)")[, 2]
keep      <- which(acc_all %in% GSM_IDS)

sample_meta <- do.call(rbind, lapply(keep, function(i) {
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
    GSM   = acc_all[i],
    title = title,
    age_weeks = as.numeric(str_extract(get_char("age"), "[0-9]+")),
    sex   = get_char("Sex"),
    stringsAsFactors = FALSE
  )
}))

# Filename prefixes on disk (matched by hand against `ls`, since the file
# prefix isn't derivable purely from the soft file's title field -- e.g.
# GSM8717008's title is "20240208-NMK12F2" but its files are prefixed
# "GSM8717008_20240208-NMK12F2U1").
file_prefix <- c(
  GSM8717006 = "GSM8717006_20210907-NMK92M-Fp1U1",
  GSM8717007 = "GSM8717007_20210908-NMK92F-Fp1U1",
  GSM8717008 = "GSM8717008_20240208-NMK12F2U1",
  GSM8717009 = "GSM8717009_20240208-NMK12M2U1"
)

sample_meta <- sample_meta %>%
  mutate(
    Age        = paste0("W", age_weeks),
    sample_id  = paste0(Age, "_", sex),
    file_prefix = file_prefix[GSM],
    matrix_dir  = DATA_DIR
  ) %>%
  arrange(age_weeks, sex) %>%
  select(GSM, sample_id, Age, age_weeks, sex, file_prefix, matrix_dir, title)

stopifnot(nrow(sample_meta) == 4)
stopifnot(!anyDuplicated(sample_meta$sample_id))

message("\n  Final metadata table:")
print(sample_meta)

out_csv <- file.path(OUT_DIR, "XeniumDev_sample_metadata.csv")
write.csv(sample_meta, out_csv, row.names = FALSE)
message("\n  Saved: ", out_csv)
