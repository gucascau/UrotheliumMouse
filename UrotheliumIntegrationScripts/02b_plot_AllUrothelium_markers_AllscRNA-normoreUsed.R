################################################################################
# 02_plot_AllUrothelium_markers.R
#
# Input : output/AllUrothelium_harmony_integrated.rds
#
# Outputs (output/AllUrothelium_markers/):
#   1. FeaturePlot  — all urothelium markers on UMAP
#   2. DotPlot      — markers × clusters (and each metadata group)
#
# Urothelium markers:
#   Krt8, Krt18, Krt19                         — pan-urothelium keratins
#   Upk1a, Upk1b, Upk2, Upk3a, Upk3b          — uroplakins
#   Krt20, Krt5, Krt14, Trp63                  — umbrella / basal markers
#   Foxa1, Gata3, Pparg                        — urothelial TFs
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
})
options(future.globals.maxSize = 8 * 1024^3)

# ── Parameters ────────────────────────────────────────────────────────────────
N_PCS        <- 50
HARMONY_DIMS <- 1:50
RESOLUTION   <- 0.5

log_mem <- function(label) {
  g  <- gc(verbose = FALSE)
  mb <- sum(g[, 2])
  message(sprintf("  [mem] %s: %.1f GB in use", label, mb / 1024))
}

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR  <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/FinalUrotheliumCells"
SCR_DIR   <- file.path(BASE_DIR, "UrotheliumScripts")
OUT_DIR   <- file.path(SCR_DIR, "output")
PLOT_DIR  <- file.path(OUT_DIR, "AllUrothelium_gated_UMAP_nosnRNA.rds")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

args     <- commandArgs(trailingOnly = TRUE)
RDS_PATH <- if (length(args) >= 1) args[[1]] else
  file.path(OUT_DIR, "AllUrothelium_harmony_integrated.rds")

# ── Marker genes ──────────────────────────────────────────────────────────────
URO_MARKERS <- c(
  "Krt8",  "Krt18", "Krt19",
  "Upk1a", "Upk1b", "Upk2", "Upk3a", "Upk3b",
  "Krt20", "Krt5",  "Krt14", "Trp63",
  "Foxa1", "Gata3", "Pparg"
)

# Metadata groups for DotPlot stratification
DOTPLOT_GROUPS <- c("seurat_clusters", "condition", "tissue", "sample_id",
                    "technology", "paper", "Categories")

# Drop conditions with too few cells or developmental origin
DROP_CONDITIONS <- c("E9To13.5Gestation", "E18_5_Kidney")

# Drop samples that should be excluded entirely (e.g. developmental datasets)
DROP_SAMPLES <- c("MultiOmicSpatialMouseKidney_SexSpecific_lifespan_Chen2025NatGenet")

# ── Load ──────────────────────────────────────────────────────────────────────
message("Loading: ", RDS_PATH)
if (!file.exists(RDS_PATH))
  stop("Input RDS not found: ", RDS_PATH)

so <- readRDS(RDS_PATH)
message(sprintf("  Loaded: %d cells × %d genes", ncol(so), nrow(so)))
log_mem("after load")

# Measure the number of cells of each metadata group (for later plotting reference), including the number of cells in each condition, tissue, sample_id, paper, technology, and seurat_clusters if they exist in the metadata.

so@meta.data %>%
  select(any_of(c("condition", "tissue", "sample_id", "paper", "technology", "seurat_clusters"))) %>%
  pivot_longer(everything(), names_to = "metadata_group", values_to = "group") %>%
  count(metadata_group, group, name = "n_cells") %>%
  group_by(metadata_group) %>%
  arrange(desc(n_cells)) %>%
  print(n = Inf)

so@meta.data$condition %>% table() %>% print()
so@meta.data$tissue %>% table() %>% print()
# echeck the NA values in tissue
so@meta.data$tissue %>% is.na() %>% table() %>% print()
# select rows of the NA values in tissue
so@meta.data %>% filter(is.na(tissue)) %>% head()
# rename the NA values in tissue to "kidney", and change the GSM3827175 and GSM3827176 samples to "bladder" based on the sample_id column
so@meta.data<- so@meta.data %>%mutate(tissue = case_when(
  grepl("GSM3827175|GSM3827176", sample_id, ignore.case = TRUE) ~ "bladder",
  is.na(tissue) ~ "kidney",
   TRUE ~ as.character(tissue)
)) 

# simpily select the sample, condition, and tissue columns, and ignore the rownames double and check the unique
so@meta.data %>% select(sample_id, condition, tissue,paper) 
df_selected <- so@meta.data %>%
  select(sample_id, condition, tissue, paper) %>%
  as.data.frame()

rownames(df_selected) <- NULL

uniquedf_selected <- df_selected %>%
  distinct()
uniquedf_selected %>% head()
uniquedf_selected %>% print()


# We noticed that some errors in the GSM3827175 and GSM3827176 samples, which are from the KUDO paper. These two samples are supposed to be from the bladder organoid, but they are labeled as "NA" in the tissue column. We will need to check the original metadata to see if we can fill in these missing values based on other columns such as sample_id or condition. If we cannot fill in the missing values, we may need to exclude these cells from our analysis or treat them as a separate category.
so@meta.data %>% filter(grepl("KUDO_GSM3827175_NMU_O_P", sample_id, ignore.case = TRUE)) %>%
head()


# check the NA values in other metadata columns
so@meta.data$tissue %>% is.na() %>% table() %>% print()
so@meta.data$sample_id %>% is.na() %>% table() %>% print()
so@meta.data$technology %>% is.na() %>% table() %>% print()
so@meta.data$condition %>% is.na() %>% table() %>% print()
so@meta.data$seurat_clusters %>% is.na() %>% table() %>% print()

so@meta.data$condition %>% table() %>% print()

so@meta.data$sample_id %>% table() %>% print()

# we create a new column called "Categories" in the metadata that seperate KidneyOrganoid, BladderOrganoid, and UreterOrganoid based on the sample_id column. We will use the following criteria:
# - If the sample_id contains "KUDO_GSM3827175_NMU_O_P" or "KUDO_GSM3827176_NMU_O_P", we will label it as "BladderOrganoid".
# - If the sample_id contains "KUDO_Yoda", KUDO_Vehicle, KUDO_GOF and KUDO_LOF , we will label it as "KidneyOrganoid"
# - If the sample_id contains "Reconstituted_ureter", we will label it as "UreterOrganoid".
# - other sample will lable as tissue paste with mouse

so@meta.data <- so@meta.data%>% mutate(Categories = case_when(
  grepl("KUDO_GSM3827175_NMU_O_P|KUDO_GSM3827176_NMU_O_P", sample_id, ignore.case = TRUE) ~ "bladderOrganoid",
  grepl("KUDO_Yoda|KUDO_Vehicle|KUDO_GOF|KUDO_LOF", sample_id, ignore.case = TRUE) ~ "kidneyOrganoid",
  grepl("MouseUreterRecon", sample_id, ignore.case = TRUE) ~ "ureterOrganoid",
  TRUE ~ paste0(tissue, "Mouse")
)) 

# add the organoid and invivo labels to the Categories column based on the sample_id column. We will use the following criteria:
so@meta.data <- so@meta.data%>% mutate(LabelClass = case_when(
  grepl("KUDO_GSM3827175_NMU_O_P|KUDO_GSM3827176_NMU_O_P", sample_id, ignore.case = TRUE) ~ "Organoid",
  grepl("KUDO_Yoda|KUDO_Vehicle|KUDO_GOF|KUDO_LOF", sample_id, ignore.case = TRUE) ~ "Organoid",
  grepl("MouseUreterRecon", sample_id, ignore.case = TRUE) ~ "Organoid",
  TRUE ~ paste0(tissue, "Invivo")
)) 

so@meta.data$Categories  %>% table()
so@meta.data$sample_id  %>% table()

# ── Fix metadata for GSE190887 (UUO_CellMeta2022, sci-RNA-seq3) ───────────────
# Cells with sample_id == "Urotherlium" originate from GSE190887 (Wu et al.,
# a mouse kidney sci-RNA-seq3 atlas covering Health, UUO, and IRI timepoints).
# Their current metadata has three errors:
#   condition  = "UUO"        (wrong — mix of Health, IRI_*, and UUO_* cells)
#   technology = "10X"        (wrong — sci-RNA-seq3; 20-mer plate-indexed barcode)
#   gsm_id     = NA
# We look up each cell's true condition via its sci-RNA-seq3 plate barcode.
GEO_META_PATH <- file.path(
  "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse",
  "UUO_CellMeta2022/GSE190887_meta_cell_type_sample.csv"
)

if (file.exists(GEO_META_PATH) &&
    any(so@meta.data$sample_id == "Urotherlium", na.rm = TRUE)) {

  message("Fixing GSE190887 ('Urotherlium') cell metadata ...")
  geo_meta <- read.csv(GEO_META_PATH, row.names = 1,
                       stringsAsFactors = FALSE)

  uro_idx  <- which(so@meta.data$sample_id == "Urotherlium")

  # Convert integrated barcode → GSE190887 barcode format
  #   "RenalUrothelium_P5001_AGTATTAGCTTCGTAGAGAA-17"
  #   → strip prefix → "P5001_AGTATTAGCTTCGTAGAGAA-17"
  #   → strip lib suffix → "P5001_AGTATTAGCTTCGTAGAGAA"
  #   → plate separator → "P5001.AGTATTAGCTTCGTAGAGAA"
  bc_raw  <- rownames(so@meta.data)[uro_idx]
  bc      <- sub("^RenalUrothelium_", "", bc_raw)
  bc      <- sub("-\\d+$",            "", bc)
  bc_geo  <- sub("^(P\\d+)_",         "\\1.", bc)

  matched_cond <- geo_meta[bc_geo, "sample"]
  n_matched    <- sum(!is.na(matched_cond))
  message(sprintf("  %d / %d cells matched to GSE190887 barcodes",
                  n_matched, length(uro_idx)))

  so@meta.data$condition[uro_idx]  <- ifelse(
    !is.na(matched_cond), matched_cond, so@meta.data$condition[uro_idx]
  )
  so@meta.data$technology[uro_idx] <- "sci-RNA-seq3"
  so@meta.data$gsm_id[uro_idx]     <- "GSE190887"

  message("  Condition breakdown after fix:")
  print(table(so@meta.data$condition[uro_idx], useNA = "ifany"))
} else {
  message("  GSE190887 metadata not found or no 'Urotherlium' cells — skipping fix.")
}
# ── Fix metadata for source_GEO (AKI, Controls) ───────────────
so@meta.data<- so@meta.data%>% mutate(condition = case_when(
  grepl("GSE209610", source_GEO, ignore.case = TRUE) ~ condition_level1,
  grepl("GSE119531", source_GEO, ignore.case = TRUE) ~ "UUO_14days",
  TRUE ~ condition
)) 

# create a new column called "FinalConditions" 
so@meta.data$condition %>% table()

# Healthy_Urothelium --> HealthyBladder
# Urothelium_Organoid --> BladderUroOganoid
# 

# -- Fix metadata for Sample_ID and orign.ident
# we have a duplicated name of Urothelium in the sample_id column, which is from the GSE190887 dataset. We will rename the sample_id of the cell to combine identity and conditions for these cells/
# We have a duplicated samples from this GSE190887 dataset and BladderNormal1
so@meta.data <- so@meta.data  %>% mutate(DuplicateSamples = case_when(
  grepl("UUO_Day10|UUO_Day14|UUO_Day2|UUO_Day4|UUO_Day6", condition) ~ "Duplicated",
  grepl("BladderNormal1|BladderNormal2", sample_id) ~ "Duplicated",
  grepl("Urotherlium", sample_id) & grepl("Health",condition) ~ "Duplicated",
  TRUE ~ "Unique"
))

# we filtered the samples that contain duplicated sample_id.
so <- subset(so, subset = DuplicateSamples == "Unique")
so@meta.data$sample_id %>% table()
so@meta.data$orig.ident %>% table()
# we change the orig.ident with the urothelum name and combine with condition
so@meta.data<- so@meta.data %>% mutate(orig.ident = case_when(
  # UUO_Day*/Healthy Urotherlium cells are already removed above; only IRI
  # Urotherlium cells remain — match all by name.
  # Note: the actual value is "Urotherlium" (with an extra 'r').
  grepl("Urotherlium", orig.ident) ~ paste0(orig.ident, "_", condition),
  grepl("D2|D1", orig.ident)  ~ paste0(orig.ident, "_", condition),
  grepl("UUO", orig.ident)  ~ "UUOKidney_14D",
  is.na(orig.ident) ~ sample_id,
  TRUE ~ orig.ident),
sample_id = case_when(
  grepl("Urotherlium", orig.ident) ~ orig.ident,
  grepl("D2|D1", orig.ident)  ~ orig.ident,
  grepl("UUO", orig.ident)  ~ "UUOKidney_14D",
  TRUE ~ sample_id)
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

so@meta.data %>% select(sample_id, orig.ident, condition, tissue, technology, paper, Categories, LabelClass,FinalConditionL2, FinalConditionL1, Finalgsm_id, Finalpaper ) %>% tail()

# check the NA values for FinalCategories
so@meta.data %>% filter(is.na(FinalConditionL1)) %>% pull(sample_id) %>% table()

so@meta.data %>% filter(sample_id %in% c("GOF","LOF","Yoda","Vehicle"))


so@meta.data %>% filter(is.na(Finalpaper)) %>% select(sample_id, orig.ident, condition, tissue, technology, paper, Categories, LabelClass,FinalConditionL2, FinalConditionL1, Finalgsm_id, Finalpaper )%>% tail()
################################################################################
# Urothelium selection criterion  (mirrors 03_extract_uro.py three-arm gate)
#
# A cell is kept when it passes ANY arm relevant to its tissue origin AND
# passes the kidney-score gate where applicable.
#
# Arm                  | All tissues | Kidney only
# ─────────────────────────────────────────────────
# Umbrella (any UPK)   |     YES     |
# Basal (>=2 Krt5/14)  |     YES     |
# Intermediate (strict)|     YES     |
# Pan-keratin (Krt8/18)|  non-kidney | NOT used — Krt18 is also in
#                      |             | kidney tubular cells, so it
#                      |             | cannot gate kidney urothelium
# KidneyEpiScore < 0.2 |             | qc kidney only (KUDO)
#
# Rationale for pan-keratin arm on bladder/ureter:
#   Ureter urothelium expresses lower uroplakins than bladder, and the strict
#   intermediate arm requires all four markers simultaneously — scRNA-seq
#   dropout removes many genuine cells.  Krt8/Krt18/Krt19 are universal
#   urothelial markers with no kidney-tubular ambiguity when no kidney cells
#   are present.
################################################################################

UPK_MARKERS    <- c("Upk1a", "Upk1b", "Upk2", "Upk3a", "Upk3b")
BASAL_MARKERS  <- c("Krt5", "Krt14", "Trp63")
BASAL_MIN_POS  <- 2L
LUMINAL_TF     <- c("Foxa1", "Gata3")
PAN_KRT        <- c("Krt8", "Krt18", "Krt19")  # non-kidney arm only

KIDNEY_EXCLUDE <- c(
  "Slc34a1", "Lrp2",    "Cubn",     # proximal tubule
  "Umod",    "Slc12a1",              # TAL
  "Slc12a3",                         # DCT
  "Aqp2",    "Aqp3",    "Scnn1g",   # collecting duct principal
  "Atp6v1b1", "Slc4a1",  "Foxi1"   # intercalated cells
)
KIDNEY_SCORE_MAX <- 0.2

# Helper: count expressed genes per cell from a gene set
n_positive <- function(so, genes) {
  genes_use <- intersect(genes, rownames(so))
  if (length(genes_use) == 0L) return(rep(0L, ncol(so)))
  expr <- GetAssayData(so, layer = "data")[genes_use, , drop = FALSE]
  Matrix::colSums(expr > 0)
}
any_positive <- function(so, genes) n_positive(so, genes) >= 1L

message("\nChecking marker availability for urothelium gate ...")
for (nm in list(UPK = UPK_MARKERS, BASAL = BASAL_MARKERS,
                LUMINAL_TF = LUMINAL_TF, KIDNEY_EXCLUDE = KIDNEY_EXCLUDE)) {
  absent <- setdiff(nm, rownames(so))
  if (length(absent) > 0)
    message(sprintf("  Missing: %s", paste(absent, collapse = ", ")))
}

# -- KidneyEpiScore via AddModuleScore ----------------------------------------
kidney_found <- intersect(KIDNEY_EXCLUDE, rownames(so))
if (length(kidney_found) == 0)
  stop("No kidney-exclusion markers found -- cannot compute KidneyEpiScore.")

message(sprintf(
  "Computing KidneyEpiScore1 from %d markers ...", length(kidney_found)
))
so <- AddModuleScore(
  so,
  features = list(kidney_found),
  name     = "KidneyEpiScore",   # Seurat appends "1" -> KidneyEpiScore1
  ctrl     = min(50L, nrow(so) - length(kidney_found)),
  seed     = 0
)
message(sprintf("  KidneyEpiScore1: min=%.3f  median=%.3f  max=%.3f",
                min(so$KidneyEpiScore1), median(so$KidneyEpiScore1),
                max(so$KidneyEpiScore1)))

# KidneyEpiScore applied ONLY to qc kidney cells (KUDO).
# scVI kidney cells were already filtered by 03_extract_uro.py; the
# AddModuleScore background shifts on a mixed-tissue dataset, inflating
# scores for cells that are genuine urothelium.
# Bladder/ureter cells have no tubular contaminants — score irrelevant.
needs_kidney_gate <- so$tissue == "kidney" & so$source == "qc"
kidney_low_mask <- (
  (!needs_kidney_gate) | (so$KidneyEpiScore1 < KIDNEY_SCORE_MAX)
)
message(sprintf(
  "  KidneyEpiScore filter: %s qc kidney cells (KUDO) only.",
  format(sum(needs_kidney_gate), big.mark = ",")
))

# -- Four-arm gate -------------------------------------------------------------
upk_mask   <- any_positive(so, UPK_MARKERS)
basal_mask <- n_positive(so, BASAL_MARKERS) >= BASAL_MIN_POS
interm_mask <-
  any_positive(so, "Epcam") &
  any_positive(so, "Krt7")  &
  any_positive(so, "Krt8")  &
  any_positive(so, LUMINAL_TF)
# Pan-keratin arm: only active for non-kidney cells (Krt18 is expressed in
# kidney proximal tubule, so it cannot safely gate kidney urothelium alone).
is_kidney <- so$tissue == "kidney"
pan_krt_mask <- (!is_kidney) & any_positive(so, PAN_KRT)

strict_uro_mask <-
  (upk_mask | basal_mask | interm_mask | pan_krt_mask) & kidney_low_mask

so$uro_umbrella_arm     <- upk_mask
so$uro_basal_arm        <- basal_mask
so$uro_intermediate_arm <- interm_mask
so$uro_pan_krt_arm      <- pan_krt_mask
so$strict_urothelium    <- strict_uro_mask

# -- Report -------------------------------------------------------------------
n_cells <- ncol(so)
report_gate <- function(label, mask) {
  n <- sum(mask)
  message(sprintf(
    "  %-52s: %s / %s (%.2f%%)",
    label,
    format(n, big.mark = ","),
    format(n_cells, big.mark = ","),
    100 * n / n_cells
  ))
}
message("\nUrothelial extraction gates (all cells):")
report_gate("any UPK > 0  (umbrella arm)",           upk_mask)
report_gate(">=2/3 Krt5/Krt14/Trp63 (basal arm)",    basal_mask)
report_gate("Epcam+Krt7+Krt8+TF (intermediate arm)", interm_mask)
report_gate("Krt8/18/19 > 0 (pan-krt, non-kidney)",  pan_krt_mask)
report_gate(
  "any arm positive",
  upk_mask | basal_mask | interm_mask | pan_krt_mask
)
report_gate(
  sprintf("KidneyEpiScore1 < %.2f (qc kidney only)", KIDNEY_SCORE_MAX),
  needs_kidney_gate & (so$KidneyEpiScore1 < KIDNEY_SCORE_MAX)
)
report_gate("FINAL: (any arm) AND kidney-score gate", strict_uro_mask)

# Per-tissue breakdown
if ("tissue" %in% colnames(so@meta.data)) {
  message("\nCells removed per tissue:")
  for (tis in sort(unique(so$tissue))) {
    tis_mask <- so$tissue == tis
    n_keep <- sum(strict_uro_mask & tis_mask)
    n_drop <- sum(!strict_uro_mask & tis_mask)
    message(sprintf(
      "  %-10s: keep %s  drop %s",
      tis,
      format(n_keep, big.mark = ","),
      format(n_drop, big.mark = ",")
    ))
  }
}

n_drop_uro <- sum(!strict_uro_mask)
if (n_drop_uro > 0) {
  message(sprintf(
    "  Removing %s non-urothelial cells; %s remain.",
    format(n_drop_uro, big.mark = ","),
    format(sum(strict_uro_mask), big.mark = ",")
  ))
  so <- so[, strict_uro_mask]
} else {
  message("  All cells passed the urothelium gate -- no cells removed.")
}

so@meta.data |>
  filter(grepl("KUDO", sample_id, ignore.case = TRUE)) |>
  group_by(sample_id) |>
  summarise(n_cells = n()) |>
  print(n = Inf)

# Drop developmental datasets if present
if ("condition" %in% colnames(so@meta.data)) {
  keep <- !(as.character(so$condition) %in% DROP_CONDITIONS)
  n_drop <- sum(!keep)
  if (n_drop > 0) {
    so <- so[, keep]
    message(sprintf("  Dropped %s developmental cells; %s cells remain.",
                    format(n_drop, big.mark = ","),
                    format(ncol(so), big.mark = ",")))
  }
}

# Drop excluded samples (developmental or otherwise unwanted)
if ("sample_id" %in% colnames(so@meta.data)) {
  keep <- !(as.character(so$sample_id) %in% DROP_SAMPLES)
  n_drop <- sum(!keep)
  if (n_drop > 0) {
    so <- so[, keep]
    message(sprintf("  Dropped %s cells from excluded samples; %s cells remain.",
                    format(n_drop, big.mark = ","),
                    format(ncol(so), big.mark = ",")))
  }
}

# ── Restore data layer for expression plots ───────────────────────────────────
# The integrated object may have had scale.data freed to save memory.
# data layer (log-norm) should still be present; if not, reload from source.
rna_layers <- Layers(so[["RNA"]])
if (!"data" %in% rna_layers) {
  message("  data layer missing — joining layers to restore ...")
  so <- JoinLayers(so)
}

# Verify markers present in the data
present <- intersect(URO_MARKERS, rownames(so))
missing <- setdiff(URO_MARKERS, rownames(so))
if (length(missing) > 0)
  message(sprintf("  Markers not found in object (skipped): %s",
                  paste(missing, collapse = ", ")))
if (length(present) == 0)
  stop("None of the requested marker genes found in the object.")

message(sprintf("  Plotting %d / %d markers: %s",
                length(present), length(URO_MARKERS),
                paste(present, collapse = ", ")))


# Save gated object for downstream use
saveRDS(so, file.path(OUT_DIR, "AllUrothelium_markers_gated.rds"))

### run the Seruat and batch removal pipeline on the gated object, and save the integrated object for downstream use. We will use the same parameters as in the original integration, but we will only use the gated cells. We will also save the integrated object with a new name to avoid overwriting the original integrated object.

################################################################################
# STEP : Finder variable features (HVGs)
################################################################################


so <- FindVariableFeatures(so, selection.method = "vst",
                               nfeatures = 3000, verbose = FALSE)
message(sprintf("  Top 10 HVGs: %s",
                paste(head(VariableFeatures(so), 10), collapse = ", ")))


################################################################################
# STEP 4: ScaleData
################################################################################

# Fill any remaining pct_mt NAs before regression
n_na <- sum(is.na(so$pct_mt))
if (n_na > 0) {
  so$pct_mt[is.na(so$pct_mt)] <- median(so$pct_mt, na.rm = TRUE)
}

message("ScaleData (regress pct_mt, HVGs only) ...")
so  <- ScaleData(so ,
                    features        = VariableFeatures(so),
                    vars.to.regress = "pct_mt",
                    verbose         = FALSE)
message("  ScaleData done")
log_mem("after ScaleData")


################################################################################
# STEP 5: RunPCA
################################################################################

message(sprintf("RunPCA (%d PCs) ...", N_PCS))
so  <- RunPCA(so , npcs = N_PCS, verbose = FALSE)
message("  PCA done")


################################################################################
# STEP 6: RunHarmony
################################################################################

# Use sample_id always; add technology if >1 level
harmony_vars <- "sample_id"
for (v in c("technology","FinalCategories")) {
  if (v %in% colnames(so@meta.data) &&
      length(unique(so@meta.data[[v]])) > 1) {
    harmony_vars <- c(harmony_vars, v)
  }
}

message(sprintf("RunHarmony (batch = %s) ...", paste(harmony_vars, collapse = " + ")))

so <- RunHarmony(
  so,
  group.by.vars    = harmony_vars,
  reduction        = "pca",
  reduction.save   = "harmony",
  plot_convergence = FALSE,
  project.dim      = FALSE,
  verbose          = FALSE,
  theta = 5
)
so[["pca"]] <- NULL
gc()
message("  Harmony done")
log_mem("after Harmony")


################################################################################
# STEP 7: FindNeighbors + FindClusters
################################################################################

message("FindNeighbors (annoy, k = 20) ...")
so <- FindNeighbors(
  so,
  reduction    = "harmony",
  dims         = HARMONY_DIMS,
  nn.method    = "annoy",
  k.param      = 20,
  annoy.metric = "euclidean",
  n.trees      = 50,
  verbose      = FALSE
)
log_mem("after FindNeighbors")

message(sprintf("FindClusters (resolution = %.2f) ...", RESOLUTION))
so <- FindClusters(so, resolution = RESOLUTION, verbose = FALSE)
message(sprintf("  Clusters: %d", length(unique(so$seurat_clusters))))
gc()


################################################################################
# STEP 8: RunUMAP
################################################################################

message("RunUMAP on Harmony embedding ...")
so <- RunUMAP(so, reduction = "harmony", dims = HARMONY_DIMS,
                  reduction.name = "umap_harmony", verbose = FALSE)
message("  UMAP done")

# save the object with the UMAP embedding
# versioned snapshot + canonical name (used by script 03 and downstream)
saveRDS(so, file.path(OUT_DIR, "AllUrothelium_postHarmony_markergated_UMAP_V1_050326.rds"))
saveRDS(so, file.path(OUT_DIR, "AllUrothelium_postHarmony_markergated_UMAP.rds"))

################################################################################
# 1. FeaturePlot — all markers on UMAP
################################################################################

message("Generating FeaturePlot ...")
n_col <- 5
n_row <- ceiling(length(present) / n_col)

fp <- FeaturePlot(
  so,
  features       = present,
  reduction      = "umap_harmony",
  ncol           = n_col,
  raster         = TRUE,
  order          = TRUE,
  cols           = c("lightgrey", "#d73027")
)

pdf(file.path(PLOT_DIR, "AllUrothelium_FeaturePlot_markers.pdf"),
    width = n_col * 4, height = n_row * 4)
print(fp)
dev.off()
message("  Saved: AllUrothelium_FeaturePlot_markers.pdf")


fp <- FeaturePlot(
  so,
  features       = present,
  reduction      = "umap_harmony",
  split.by = "Categories",
  ncol           = n_col,
  raster         = TRUE,
  order          = TRUE,
  cols           = c("lightgrey", "#d73027")
)

pdf(file.path(PLOT_DIR, "AllUrothelium_FeaturePlot_markers_splitCategories.pdf"),
    width = n_col * 8, height = n_row * 4)
print(fp)
dev.off()
message("  Saved: AllUrothelium_FeaturePlot_markers_splitCategories.pdf")




################################################################################
# 2. Individual FeaturePlot PDFs (one per gene, cleaner resolution)
################################################################################

message("Generating individual gene FeaturePlots ...")
for (gene in present) {
  p <- FeaturePlot(
    so,
    features  = gene,
    reduction = "umap_harmony",
    raster    = TRUE,
    order     = TRUE,
    cols      = c("lightgrey", "#d73027")
  ) + ggtitle(gene)
  ggsave(
    file.path(PLOT_DIR, sprintf("FeaturePlot_%s.pdf", gene)),
    plot = p, width = 6, height = 5
  )
}
message(sprintf("  Saved %d individual FeaturePlot PDFs", length(present)))


################################################################################
# 3. DotPlot — markers × each metadata grouping
################################################################################

message("Generating DotPlots ...")

for (grp in DOTPLOT_GROUPS) {
  if (!grp %in% colnames(so@meta.data)) next

  # Drop NA / empty-string levels and levels with < 3 cells to avoid
  # a Seurat DotPlot crash (0-row data.frame) on sparse groups.
  grp_vals <- as.character(so@meta.data[[grp]])
  grp_vals[is.na(grp_vals) | grp_vals == ""] <- NA
  valid_levels <- names(which(table(grp_vals) >= 3))
  if (length(valid_levels) == 0) {
    message(sprintf("  Skipping DotPlot for %s: no groups with >= 3 cells", grp))
    next
  }
  so_sub <- so[, !is.na(grp_vals) & grp_vals %in% valid_levels]
  so_sub@meta.data[[grp]] <- droplevels(factor(so_sub@meta.data[[grp]],
                                               levels = valid_levels))
  n_dropped <- ncol(so) - ncol(so_sub)
  if (n_dropped > 0)
    message(sprintf("  DotPlot %s: dropped %d cells in sparse/NA groups", grp, n_dropped))

  n_levels <- length(valid_levels)

  # Scale figure height to number of groups
  height <- max(4, 2 + n_levels * 0.35)
  width  <- max(8, 2 + length(present) * 0.6)

  p <- DotPlot(
    so_sub,
    features  = present,
    group.by  = grp,
    cols      = c("lightgrey", "#08519c"),
    dot.scale = 6
  ) +
    RotatedAxis() +
    labs(title = sprintf("Urothelium markers - grouped by %s", grp)) +
    theme(axis.text.x = element_text(size = 9),
          axis.text.y = element_text(size = 9))

  fname <- sprintf("AllUrothelium_DotPlot_by_%s.pdf", grp)
  ggsave(file.path(PLOT_DIR, fname), plot = p,
         width = width, height = height)
  message(sprintf("  Saved: %s", fname))
}


################################################################################
# 4. Combined overview: DimPlot clusters + DotPlot by cluster side-by-side
################################################################################

message("Generating combined overview plot ...")

p_dim <- DimPlot(so, group.by = "seurat_clusters", reduction = "umap_harmony",
                 label = TRUE, repel = TRUE, raster = TRUE) +
  ggtitle("Clusters") + NoLegend()

p_dot <- DotPlot(so, features = present, group.by = "seurat_clusters",
                 cols = c("lightgrey", "#08519c"), dot.scale = 5) +
  RotatedAxis() +
  labs(title = "Markers by cluster") +
  theme(axis.text.x = element_text(size = 8))

pdf(file.path(PLOT_DIR, "AllUrothelium_overview_clusters_dotplot.pdf"),
    width = 22, height = 10)
print(p_dim | p_dot)
dev.off()
message("  Saved: AllUrothelium_overview_clusters_dotplot.pdf")

# DotPlot split by tissue (side-by-side, one panel per tissue)
if ("tissue" %in% colnames(so@meta.data)) {
  tissues <- sort(unique(so@meta.data$tissue))
  dot_by_tissue <- lapply(tissues, function(tis) {
    so_sub <- so[, so$tissue == tis]
    if (ncol(so_sub) < 10) return(NULL)
    DotPlot(so_sub, features = present, group.by = "seurat_clusters",
            cols = c("lightgrey", "#08519c"), dot.scale = 5) +
      RotatedAxis() +
      ggtitle(sprintf("Tissue: %s", tis)) +
      theme(axis.text.x = element_text(size = 7))
  })
  dot_by_tissue <- Filter(Negate(is.null), dot_by_tissue)

  if (length(dot_by_tissue) > 0) {
    pdf(file.path(PLOT_DIR, "AllUrothelium_DotPlot_byTissue_panels.pdf"),
        width = 20, height = 8 * length(dot_by_tissue))
    print(wrap_plots(dot_by_tissue, ncol = 1))
    dev.off()
    message("  Saved: AllUrothelium_DotPlot_byTissue_panels.pdf")
  }
}

message("Marker plots complete. Output: ", PLOT_DIR)



################################################################################
# 4. Combined overview: DimPlot plots by each categories
################################################################################
message("Generating Dimplots ...")

for (grp in DOTPLOT_GROUPS) {
  if (!grp %in% colnames(so@meta.data)) next

  n_levels <- length(unique(so@meta.data[[grp]]))

  # Scale figure height to number of groups
  height <- max(4, 2 + n_levels * 0.35)
  width  <- max(8, 2 + length(present) * 0.6)
  p <- DimPlot(
    so, group.by = grp, reduction = "umap_harmony",
    label = FALSE, repel = TRUE, raster = TRUE
  ) +
    NoLegend() +
    labs(title = sprintf("Umap — grouped by %s", grp)) +
    theme(
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9)
    )

  fname <- sprintf("AllUrothelium_DimPlot_by_%s.pdf", grp)
  ggsave(file.path(PLOT_DIR, fname), plot = p,
         width = width, height = width)
  message(sprintf("  Saved: %s", fname))
}



for (grp in DOTPLOT_GROUPS) {
  if (!grp %in% colnames(so@meta.data)) next

  n_levels <- length(unique(so@meta.data[[grp]]))

  # Scale figure height to number of groups
  height <- max(4, 2 + n_levels * 0.35)
  width  <- max(15, 2 + length(present) * 0.6)
  p <- DimPlot(
    so, group.by = "seurat_clusters", split.by = grp,
    label = FALSE, repel = TRUE, raster = TRUE
  ) +
    NoLegend() +
    labs(title = sprintf("Umap — grouped by %s", grp)) +
    theme(
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9)
    )

  fname <- sprintf("AllUrothelium_DimPlot_Splitby_%s.pdf", grp)
  ggsave(file.path(PLOT_DIR, fname), plot = p,
         width = width, height = height)
  message(sprintf("  Saved: %s", fname))
}




