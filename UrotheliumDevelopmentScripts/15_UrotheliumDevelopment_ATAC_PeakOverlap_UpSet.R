################################################################################
# 15_UrotheliumDevelopment_ATAC_PeakOverlap_UpSet.R
#
# Figure J. Peak overlap across developmental stages -- an UpSet plot showing
#   how many peaks are stage-specific vs. shared across two or more stages,
#   plus the summary CSV of every combination's peak count.
#
# Peak "presence" in a stage reuses the exact counts-matrix proxy introduced
# in 14_UrotheliumDevelopment_ATAC_QC_Figures.R (Figure I): a peak counts as
# present in a stage if it has nonzero accessibility (count > 0) in at least
# MIN_CELLS_DETECTED cells of that stage. Same caveat as script 14 -- the
# original MACS peak_called_in field (CallPeaks(group.by="Age") in script 08)
# did not survive into UrotheliumOnly_RNA_ATAC.rds (FeatureMatrix() /
# CreateChromatinAssay() rebuild the assay from a plain counts matrix with no
# mcols), so this is a coverage-based proxy, not the original per-stage MACS
# call. Confirmed by direct user choice when this same issue came up for
# Figure I.
#
# ComplexHeatmap::make_comb_mat(..., mode = "distinct") is used rather than
# "union"/"intersect" -- "distinct" assigns each peak to exactly one
# combination (the exact set of stages it's present in), which is the
# classic UpSet semantics and the only mode where "specific" (single-stage)
# vs. "overlapping" (multi-stage) bars partition the peak set without
# double-counting.
#
# A small number of peaks in the merged set clear no stage's >= 3-cell
# threshold at all (the "all-FALSE" combination) -- these are excluded from
# the plotted matrix (kept only in the summary CSV) since they aren't
# actually assignable to any stage and would otherwise show up as a
# meaningless "0000...0" bar.
#
# Input:  UrotheliumDevelopmentScripts/output/UrotheliumOnly_RNA_ATAC.rds
# Output: FigJ_UrotheliumDevelopment_ATAC_PeakOverlap_UpSet.pdf
#         UrotheliumDevelopment_ATAC_PeakOverlap_Combinations.csv
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(Matrix)
  library(ComplexHeatmap)
  library(dplyr)
})

SCRIPT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/UrotheliumDevelopmentScripts"
OUT_DIR    <- file.path(SCRIPT_DIR, "output")

STAGE_ORDER <- c("E16.5", "P0", "W3", "W12", "W52", "W92")
MIN_CELLS_DETECTED <- 3  # same threshold as script 14's Figure I

# ── Load the combined RNA+ATAC Uro object ───────────────────────────────────
message("==> Loading Uro RNA+ATAC object ...")
in_rds <- file.path(OUT_DIR, "UrotheliumOnly_RNA_ATAC.rds")
if (!file.exists(in_rds)) {
  stop("Missing ", in_rds, " -- run 08_UrotheliumDevelopment_ATAC_Build.R first.")
}
uro <- readRDS(in_rds)
uro$Age <- factor(uro$Age, levels = STAGE_ORDER)
message(sprintf("  %d cells, %d ATAC peaks", ncol(uro), nrow(uro[["ATAC"]])))

# ── Per-stage peak presence matrix (peaks x stages, logical) ───────────────
message(sprintf("\n==> Building peak presence matrix (count > 0 in >= %d cells/stage) ...",
                 MIN_CELLS_DETECTED))
counts_mat <- GetAssayData(uro, assay = "ATAC", layer = "counts")

presence_mat <- sapply(STAGE_ORDER, function(st) {
  cells_st <- colnames(uro)[uro$Age == st]
  Matrix::rowSums(counts_mat[, cells_st, drop = FALSE] > 0) >= MIN_CELLS_DETECTED
})
rownames(presence_mat) <- rownames(counts_mat)
message("  Peaks present per stage:")
print(colSums(presence_mat))

n_undetected <- sum(rowSums(presence_mat) == 0)
message(sprintf("  Peaks present in NO stage at this threshold: %d / %d (excluded from the plot)",
                 n_undetected, nrow(presence_mat)))

# ── Combination matrix + summary CSV (all combinations, including 0-stage) ─
comb_mat_full <- make_comb_mat(presence_mat, mode = "distinct")
comb_summary <- data.frame(
  combination = comb_name(comb_mat_full),
  degree      = comb_degree(comb_mat_full),
  n_peaks     = comb_size(comb_mat_full)
) %>% arrange(desc(n_peaks))
message("\n  Combination summary (top 10 by size):")
print(head(comb_summary, 10))
write.csv(comb_summary, file.path(OUT_DIR, "UrotheliumDevelopment_ATAC_PeakOverlap_Combinations.csv"),
          row.names = FALSE)

# ── UpSet plot (0-stage / all-FALSE combination excluded) ──────────────────
message("\n==> Building Figure J (UpSet plot) ...")
comb_mat_plot <- comb_mat_full[comb_degree(comb_mat_full) > 0]

upset_plot <- UpSet(
  comb_mat_plot,
  comb_order       = order(-comb_size(comb_mat_plot)),
  set_order        = STAGE_ORDER,
  top_annotation   = upset_top_annotation(comb_mat_plot, add_numbers = TRUE),
  right_annotation = upset_right_annotation(comb_mat_plot, add_numbers = TRUE),
  row_names_gp     = grid::gpar(fontsize = 10)
)

pdf(file.path(OUT_DIR, "FigJ_UrotheliumDevelopment_ATAC_PeakOverlap_UpSet.pdf"),
    width = 12, height = 6)
ht <- draw(upset_plot,
           column_title = sprintf(
             "Figure J. snATAC-seq peak overlap across development (n = %d peaks scored, %d excluded as detected in no stage)",
             sum(comb_size(comb_mat_plot)), n_undetected))
dev.off()
message("  Saved: FigJ_UrotheliumDevelopment_ATAC_PeakOverlap_UpSet.pdf")

message("\n==> Done.")
