################################################################################
# calc_uro_proportion_organoids.R
#
# Compute urothelium cell proportions for organoid datasets and combine with
# kidney/bladder results from calc_uro_proportion_h5ad.py.
#
# Strategy:
#   Small organoid _qc.rds (< 1 GB):
#     - Load RDS, count total cells, apply urothelium marker filter.
#   KudoUUOUrothelium_qc.rds (46 GB):
#     - Too large for routine marker-based filtering on standard nodes.
#     - Total cells inferred from AllUrothelium_harmony_integrated.rds (KUDO_*
#       samples) because all Kudo _qc cells enter integration unfiltered.
#     - Proportion = 100 % by pipeline design.
#
# Outputs:
#   output/uro_proportion/uro_proportion_organoids.csv
#   output/uro_proportion/uro_proportion_combined.csv
#   output/uro_proportion/uro_proportion_plots.pdf
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
})

options(future.globals.maxSize = 16 * 1024^3)

log_mem <- function(label) {
  g  <- gc(verbose = FALSE)
  mb <- sum(g[, 2])
  message(sprintf("  [mem] %s: %.1f GB", label, mb / 1024))
}

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR     <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse"
FINAL_DIR    <- file.path(BASE_DIR, "FinalUrotheliumCells")
INTEG_RDS    <- file.path(FINAL_DIR, "UrotheliumScripts/output",
                           "AllUrothelium_harmony_integrated.rds")

OUT_DIR <- file.path(BASE_DIR,
  "UsedSingleCells/UrotheliumIntegrationScripts/output/uro_proportion")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Small organoid files (< 1 GB — applied with marker filter)
SMALL_ORGANOID_RDS <- list(
  list(
    path       = file.path(FINAL_DIR, "BladderHomogenate1_qc.rds"),
    sample_id  = "BladderHomogenate1",
    tissue     = "bladder",
    condition  = "BladderOrganoid",
    technology = "10X_v2",
    paper      = "PMID_31462402",
    gsm_id     = "GSM3723360",
    data_type  = "organoid"
  ),
  list(
    path       = file.path(FINAL_DIR, "BladderHomogenate2_qc.rds"),
    sample_id  = "BladderHomogenate2",
    tissue     = "bladder",
    condition  = "BladderOrganoid",
    technology = "10X_v2",
    paper      = "PMID_31462402",
    gsm_id     = "GSM3723361",
    data_type  = "organoid"
  ),
  list(
    path       = file.path(FINAL_DIR, "MouseUreterRecon1_qc.rds"),
    sample_id  = "MouseUreterRecon1",
    tissue     = "ureter",
    condition  = "UreterOrganoid",
    technology = "10X",
    paper      = "PMID_40541956",
    gsm_id     = "GSM8635363",
    data_type  = "organoid"
  )
)

# ── Urothelium marker filter ───────────────────────────────────────────────────
# Three-arm OR gate (mirrors 06_extract_uro.py / 03_extract_uro.py)
#   umbrella  : any Upk1a/Upk1b/Upk2/Upk3a/Upk3b > 0
#   basal     : ≥2 of Krt5/Krt14/Trp63 > 0
#   intermediate: Epcam > 0 AND Krt7 > 0 AND Krt8 > 0 AND any of Foxa1/Gata3 > 0

UPK_MARKERS      <- c("Upk1a", "Upk1b", "Upk2", "Upk3a", "Upk3b")
BASAL_MARKERS    <- c("Krt5",  "Krt14", "Trp63")
LUMINAL_TF       <- c("Foxa1", "Gata3")
INTERMEDIATE_REQ <- c("Epcam", "Krt7",  "Krt8")
BASAL_MIN        <- 2L

is_expressed <- function(counts_mat, genes) {
  present <- intersect(genes, rownames(counts_mat))
  if (length(present) == 0) return(rep(FALSE, ncol(counts_mat)))
  Matrix::colSums(counts_mat[present, , drop = FALSE] > 0) > 0
}

n_expressed <- function(counts_mat, genes) {
  present <- intersect(genes, rownames(counts_mat))
  if (length(present) == 0) return(rep(0L, ncol(counts_mat)))
  as.integer(Matrix::colSums(counts_mat[present, , drop = FALSE] > 0))
}

apply_uro_filter <- function(so) {
  DefaultAssay(so) <- "RNA"
  cnt <- GetAssayData(so, layer = "counts")

  umbrella     <- is_expressed(cnt, UPK_MARKERS)
  basal        <- n_expressed(cnt, BASAL_MARKERS) >= BASAL_MIN
  inter_core   <- is_expressed(cnt, INTERMEDIATE_REQ)
  inter_req    <- {
    present_all <- intersect(INTERMEDIATE_REQ, rownames(cnt))
    if (length(present_all) < length(INTERMEDIATE_REQ)) {
      message(sprintf("    Missing intermediate markers: %s",
                      paste(setdiff(INTERMEDIATE_REQ, rownames(cnt)), collapse = ", ")))
    }
    Matrix::colSums(
      cnt[intersect(INTERMEDIATE_REQ, rownames(cnt)), , drop = FALSE] > 0
    ) == length(intersect(INTERMEDIATE_REQ, rownames(cnt)))
  }
  inter_tf   <- is_expressed(cnt, LUMINAL_TF)
  intermediate <- inter_req & inter_tf

  umbrella | basal | intermediate
}

# ── Process small organoid files ──────────────────────────────────────────────
message("\n=== Processing small organoid _qc.rds files ===")

organoid_rows <- lapply(SMALL_ORGANOID_RDS, function(entry) {
  sid <- entry$sample_id
  message(sprintf("\n  [%s] Loading %s ...", sid, basename(entry$path)))
  if (!file.exists(entry$path)) {
    warning("File not found: ", entry$path)
    return(NULL)
  }
  so <- readRDS(entry$path)
  total <- ncol(so)
  message(sprintf("    %d total cells", total))
  log_mem(sprintf("after loading %s", sid))

  uro_mask <- tryCatch(
    apply_uro_filter(so),
    error = function(e) {
      message("    Marker filter failed: ", conditionMessage(e))
      rep(NA, total)
    }
  )
  rm(so); gc(verbose = FALSE)

  uro_cnt <- if (anyNA(uro_mask)) NA_integer_ else sum(uro_mask)
  uro_pct <- if (is.na(uro_cnt)) NA_real_ else round(100 * uro_cnt / total, 4)

  message(sprintf("    Urothelium cells: %s / %d  (%s%%)",
                  if (is.na(uro_cnt)) "NA" else format(uro_cnt, big.mark = ","),
                  total,
                  if (is.na(uro_pct)) "NA" else uro_pct))

  data.frame(
    dataset_type  = "organoid",
    sample_id     = sid,
    tissue        = entry$tissue,
    condition     = entry$condition,
    technology    = entry$technology,
    paper         = entry$paper,
    gsm_id        = entry$gsm_id,
    total_cells   = total,
    uro_cells     = uro_cnt,
    uro_proportion = if (is.na(uro_cnt)) NA_real_ else round(uro_cnt / total, 6),
    uro_pct       = uro_pct,
    note          = "marker filter applied",
    stringsAsFactors = FALSE
  )
})
organoid_df <- dplyr::bind_rows(Filter(Negate(is.null), organoid_rows))


# ── KudoUUOUrothelium: infer counts from AllUrothelium integration ─────────────
message("\n=== KudoUUOUrothelium: inferring counts from AllUrothelium_harmony_integrated.rds ===")

if (!file.exists(INTEG_RDS)) stop("Integration RDS not found: ", INTEG_RDS)

message("  Loading AllUrothelium_harmony_integrated.rds (1.7 GB) ...")
integ <- readRDS(INTEG_RDS)
message(sprintf("  Loaded: %d cells total", ncol(integ)))
log_mem("after loading integrated object")

# All KUDO_* samples — these are from KudoUUOUrothelium_qc.rds (no marker filter)
kudo_mask <- grepl("^KUDO_", integ@meta.data$sample_id)
kudo_meta <- integ@meta.data[kudo_mask, ]
message(sprintf("  KUDO cells in integration: %d", sum(kudo_mask)))

if (sum(kudo_mask) > 0) {
  kudo_counts <- kudo_meta %>%
    group_by(sample_id) %>%
    summarise(
      tissue      = first(tissue),
      condition   = first(condition),
      technology  = first(technology),
      paper       = ifelse("paper" %in% names(kudo_meta), first(paper), NA_character_),
      gsm_id      = ifelse("gsm_id" %in% names(kudo_meta), first(gsm_id), NA_character_),
      uro_cells   = n(),
      .groups     = "drop"
    ) %>%
    mutate(
      dataset_type  = "organoid",
      total_cells   = uro_cells,
      uro_proportion = 1.0,
      uro_pct       = 100.0,
      note          = paste0("counts from AllUrothelium_harmony_integrated.rds; ",
                             "all Kudo _qc cells enter integration unfiltered; ",
                             "marker filter not applied (46 GB file requires high-mem node)")
    )

  # Per-sub-sample rows
  kudo_rows <- kudo_counts
  message("\n  KUDO cell counts by sample_id:")
  for (i in seq_len(nrow(kudo_rows))) {
    message(sprintf("    %-35s : %d cells", kudo_rows$sample_id[i], kudo_rows$uro_cells[i]))
  }

  # Also add an aggregate row for the parent KudoUUOUrothelium
  kudo_total_row <- data.frame(
    dataset_type  = "organoid",
    sample_id     = "KudoUUOUrothelium (aggregate)",
    tissue        = "kidney",
    condition     = "KidneyOrganoid+scRNAseq",
    technology    = "PIPseq+sci-RNA-seq3",
    paper         = "Kudo2026+PMID_36265491",
    gsm_id        = "KUDORDS",
    total_cells   = sum(kudo_rows$uro_cells),
    uro_cells     = sum(kudo_rows$uro_cells),
    uro_proportion = 1.0,
    uro_pct       = 100.0,
    note          = kudo_rows$note[1],
    stringsAsFactors = FALSE
  )

  organoid_df <- bind_rows(organoid_df, kudo_rows, kudo_total_row)
}
rm(integ); gc(verbose = FALSE)


# ── Write organoid CSV ────────────────────────────────────────────────────────
organoid_out <- file.path(OUT_DIR, "uro_proportion_organoids.csv")
write.csv(organoid_df, organoid_out, row.names = FALSE)
message(sprintf("\n  Saved: %s", organoid_out))
print(organoid_df)


# ── Combine kidney + bladder + organoid ───────────────────────────────────────
message("\n=== Combining all results ===")

kidney_csv  <- file.path(OUT_DIR, "uro_proportion_kidney.csv")
bladder_csv <- file.path(OUT_DIR, "uro_proportion_bladder.csv")

missing <- c(
  if (!file.exists(kidney_csv))  kidney_csv,
  if (!file.exists(bladder_csv)) bladder_csv
)
if (length(missing) > 0) {
  warning("The following Python output CSVs were not found — run ",
          "calc_uro_proportion_h5ad.py first:\n  ",
          paste(missing, collapse = "\n  "))
}

all_dfs <- list()
if (file.exists(kidney_csv)) {
  k <- read.csv(kidney_csv, stringsAsFactors = FALSE)
  k$note <- "marker filter applied (kidney epi score)"
  all_dfs[["kidney"]] <- k
}
if (file.exists(bladder_csv)) {
  b <- read.csv(bladder_csv, stringsAsFactors = FALSE)
  b$note <- "marker filter applied"
  all_dfs[["bladder"]] <- b
}
all_dfs[["organoid"]] <- organoid_df

# Standardise columns and bind
std_cols <- c("dataset_type", "sample_id", "tissue", "condition",
              "technology", "paper", "gsm_id",
              "total_cells", "uro_cells", "uro_proportion", "uro_pct", "note")

bind_standard <- function(df) {
  for (col in std_cols) {
    if (!col %in% names(df)) df[[col]] <- NA
  }
  df[, std_cols]
}

combined <- bind_rows(lapply(all_dfs, bind_standard))
combined$uro_pct <- round(combined$uro_proportion * 100, 4)

combined_out <- file.path(OUT_DIR, "uro_proportion_combined.csv")
write.csv(combined, combined_out, row.names = FALSE)
message(sprintf("\n  Saved: %s", combined_out))

message("\n  Summary:")
print(
  combined %>%
    group_by(dataset_type, tissue) %>%
    summarise(
      n_samples   = n(),
      total_cells = sum(total_cells, na.rm = TRUE),
      uro_cells   = sum(uro_cells,   na.rm = TRUE),
      mean_pct    = round(mean(uro_pct, na.rm = TRUE), 2),
      .groups     = "drop"
    )
)


# ── Plots ─────────────────────────────────────────────────────────────────────
message("\n=== Generating plots ===")

plot_data <- combined %>%
  filter(!grepl("aggregate", sample_id)) %>%
  mutate(
    dataset_label  = paste0(sample_id, "\n(", dataset_type, ")"),
    uro_pct_plot   = ifelse(is.na(uro_pct), 100, uro_pct),
    tissue_f       = factor(tissue, levels = c("kidney", "bladder", "ureter"))
  )

# 1. Bar chart of uro_pct per sample (faceted by dataset_type)
p1 <- ggplot(plot_data, aes(x = reorder(sample_id, uro_pct_plot),
                             y = uro_pct_plot,
                             fill = tissue_f)) +
  geom_col(width = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", uro_pct_plot)),
            hjust = -0.1, size = 2.2) +
  coord_flip() +
  facet_wrap(~ dataset_type, scales = "free_y", ncol = 1) +
  scale_fill_manual(
    values = c(kidney = "#4393C3", bladder = "#D6604D", ureter = "#74C476"),
    na.value = "grey70"
  ) +
  labs(
    title  = "Urothelium proportion per dataset",
    x      = NULL,
    y      = "Urothelium cells (%)",
    fill   = "Tissue",
    caption = paste0(
      "Kidney/bladder: proportion derived from marker-based extraction ",
      "(umbrella / basal / intermediate arm OR gate).\n",
      "Organoid: marker filter applied where feasible; ",
      "KudoUUOUrothelium reported as 100% (used directly in integration)."
    )
  ) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom",
        plot.caption = element_text(size = 7, hjust = 0))

# 2. Stacked bar: total vs uro cells per dataset
plot_stack <- combined %>%
  filter(!grepl("aggregate", sample_id), !is.na(total_cells)) %>%
  mutate(
    non_uro_cells = total_cells - uro_cells,
    tissue_f = factor(tissue, levels = c("kidney", "bladder", "ureter"))
  ) %>%
  select(sample_id, dataset_type, tissue_f, uro_cells, non_uro_cells) %>%
  pivot_longer(cols = c(uro_cells, non_uro_cells),
               names_to = "cell_type", values_to = "count") %>%
  mutate(cell_type = factor(cell_type,
                             levels = c("non_uro_cells", "uro_cells"),
                             labels = c("Non-urothelium", "Urothelium")))

p2 <- ggplot(plot_stack,
             aes(x = reorder(sample_id, count),
                 y = count, fill = cell_type)) +
  geom_col(width = 0.8) +
  coord_flip() +
  facet_wrap(~ dataset_type, scales = "free", ncol = 1) +
  scale_fill_manual(values = c("Urothelium" = "#2166AC", "Non-urothelium" = "#D1E5F0")) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Cell counts per dataset",
    x     = NULL, y = "Cell count", fill = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom")

# 3. Summary donut / bar by tissue
tissue_summary <- combined %>%
  filter(!grepl("aggregate", sample_id), !is.na(total_cells)) %>%
  group_by(tissue) %>%
  summarise(
    total_cells = sum(total_cells, na.rm = TRUE),
    uro_cells   = sum(uro_cells,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(uro_pct = round(100 * uro_cells / total_cells, 2))

p3 <- ggplot(tissue_summary,
             aes(x = reorder(tissue, uro_pct), y = uro_pct,
                 fill = tissue)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%\n(%s uro / %s total)",
                                uro_pct,
                                scales::comma(uro_cells),
                                scales::comma(total_cells))),
            vjust = -0.3, size = 3) +
  scale_fill_manual(
    values = c(kidney = "#4393C3", bladder = "#D6604D", ureter = "#74C476"),
    guide = "none"
  ) +
  ylim(0, max(tissue_summary$uro_pct) * 1.25) +
  labs(
    title = "Aggregate urothelium proportion by tissue",
    x = NULL, y = "Urothelium cells (%)"
  ) +
  theme_bw(base_size = 10)

# Save PDF
pdf_out <- file.path(OUT_DIR, "uro_proportion_plots.pdf")
tryCatch({
  pdf(pdf_out, width = 14, height = max(10, nrow(plot_data) * 0.35 + 4))
  print(p1)
  print(p2)
  print(p3)
  dev.off()
  message(sprintf("  Saved: %s", pdf_out))
}, error = function(e) {
  try(dev.off(), silent = TRUE)
  warning("PDF failed: ", conditionMessage(e))
})

message("\n=== calc_uro_proportion_organoids.R complete ===")
message(sprintf("  Combined CSV : %s", combined_out))
message(sprintf("  Plots PDF    : %s", pdf_out))
