#!/usr/bin/env Rscript
# 09_deseq2_deg.R — Pseudobulk DEG analysis with DESeq2.
#
# Reads the pseudobulk count matrix and metadata produced by 08_pseudobulk.py,
# then runs DESeq2 for each cell type comparing the specified conditions.
#
# Outputs (in output/pseudobulk/DEG/)
# ------------------------------------
#   <cell_type>_<contrast>.csv  — full DESeq2 results (all genes, ranked by padj)
#   volcano_<cell_type>_<contrast>.pdf  — volcano plot
#   summary.csv  — number of sig DEGs per cell type × contrast
#
# Usage (edit CONTRASTS and parameters below, then run via SLURM or interactively)
# ---------------------------------------------------------------------------------
#   Rscript 09_deseq2_deg.R

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
})

# ── Parameters ─────────────────────────────────────────────────────────────────
BASE_DIR  <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
PB_DIR    <- file.path(BASE_DIR, "RenalUrotheliumScripts", "output", "pseudobulk")
DEG_DIR   <- file.path(PB_DIR, "DEG")
dir.create(DEG_DIR, showWarnings = FALSE, recursive = TRUE)

# Reference condition (denominator in fold-change)
REF_CONDITION <- "Healthy"

# Comparisons: each entry is c(numerator, denominator)
# Change these to match your experimental design
CONTRASTS <- list(
  c("UUO",     "Healthy"),
  c("rUUO",    "Healthy"),
  c("TET2UUO", "Healthy"),
  c("UUO",     "rUUO")
)

# Minimum pseudobulk replicates required per group in a contrast
MIN_REPLICATES <- 3

# Significance thresholds
PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1.0     # |log2FC| ≥ 1 (2-fold)

# ── Load pseudobulk data ───────────────────────────────────────────────────────
cat("Loading pseudobulk counts and metadata ...\n")
counts <- read.csv(file.path(PB_DIR, "counts.csv"),
                   row.names = 1, check.names = FALSE)
meta   <- read.csv(file.path(PB_DIR, "meta.csv"),
                   row.names = 1, check.names = FALSE)

# Align columns
stopifnot(all(colnames(counts) == rownames(meta)))
cat(sprintf("  %d genes × %d pseudobulk samples\n", nrow(counts), ncol(counts)))
cat("  Conditions:", paste(unique(meta$condition), collapse = ", "), "\n")
cat("  Cell types:", paste(unique(meta$cell_type),  collapse = ", "), "\n\n")

# ── Helper: run DESeq2 for one cell type × one contrast ───────────────────────
run_deseq2 <- function(ct, cond_num, cond_den) {
  tag <- paste0(ct, "__", cond_num, "_vs_", cond_den)

  # Subset to this cell type and relevant conditions
  sel  <- meta$cell_type == ct & meta$condition %in% c(cond_num, cond_den)
  meta_sub   <- meta[sel, , drop = FALSE]
  counts_sub <- counts[, rownames(meta_sub), drop = FALSE]

  n_num <- sum(meta_sub$condition == cond_num)
  n_den <- sum(meta_sub$condition == cond_den)

  if (n_num < MIN_REPLICATES || n_den < MIN_REPLICATES) {
    cat(sprintf("  [%s] SKIP — %s n=%d, %s n=%d (need ≥%d each)\n",
                tag, cond_num, n_num, cond_den, n_den, MIN_REPLICATES))
    return(NULL)
  }

  cat(sprintf("  [%s] %s(n=%d) vs %s(n=%d) ...\n",
              ct, cond_num, n_num, cond_den, n_den))

  # Filter lowly expressed genes: keep genes with ≥ 10 counts in ≥ 2 samples
  keep        <- rowSums(counts_sub >= 10) >= 2
  counts_filt <- counts_sub[keep, , drop = FALSE]
  cat(sprintf("    Genes after low-count filter: %d / %d\n",
              sum(keep), nrow(counts_sub)))

  meta_sub$condition <- factor(meta_sub$condition,
                               levels = c(cond_den, cond_num))  # den = reference

  dds <- tryCatch(
    DESeqDataSetFromMatrix(
      countData = as.matrix(counts_filt),
      colData   = meta_sub,
      design    = ~ condition
    ),
    error = function(e) { cat("    DESeqDataSetFromMatrix failed:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(dds)) return(NULL)

  dds <- tryCatch(
    DESeq(dds, quiet = TRUE),
    error = function(e) { cat("    DESeq() failed:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(dds)) return(NULL)

  res <- results(dds,
                 contrast  = c("condition", cond_num, cond_den),
                 alpha     = PADJ_CUTOFF)
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    arrange(padj, desc(abs(log2FoldChange)))

  # Save full results
  out_csv <- file.path(DEG_DIR, paste0(tag, ".csv"))
  write.csv(res_df, out_csv, row.names = FALSE)

  n_up   <- sum(res_df$padj < PADJ_CUTOFF & res_df$log2FoldChange >  LFC_CUTOFF, na.rm = TRUE)
  n_down <- sum(res_df$padj < PADJ_CUTOFF & res_df$log2FoldChange < -LFC_CUTOFF, na.rm = TRUE)
  cat(sprintf("    Sig up: %d  |  Sig down: %d  (padj<%.2f, |LFC|>%.1f)\n",
              n_up, n_down, PADJ_CUTOFF, LFC_CUTOFF))

  # ── Volcano plot ─────────────────────────────────────────────────────────────
  plot_df <- res_df %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(
      sig      = padj < PADJ_CUTOFF & abs(log2FoldChange) >= LFC_CUTOFF,
      direction = case_when(
        sig & log2FoldChange >  LFC_CUTOFF ~ "Up",
        sig & log2FoldChange < -LFC_CUTOFF ~ "Down",
        TRUE                               ~ "NS"
      ),
      label = ifelse(sig & abs(log2FoldChange) >= 2 &
                       -log10(padj) >= quantile(-log10(padj[sig]), 0.8, na.rm = TRUE),
                     gene, NA)
    )

  gg <- ggplot(plot_df, aes(log2FoldChange, -log10(padj), colour = direction)) +
    geom_point(size = 0.6, alpha = 0.6) +
    scale_colour_manual(values = c(Up = "#d62728", Down = "#1f77b4", NS = "grey70")) +
    geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed",
               colour = "grey40", linewidth = 0.4) +
    geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = "dashed",
               colour = "grey40", linewidth = 0.4) +
    ggrepel::geom_text_repel(aes(label = label), size = 2.5,
                              max.overlaps = 20, na.rm = TRUE) +
    labs(title   = paste0(ct, "  |  ", cond_num, " vs ", cond_den),
         x       = expression(log[2]~fold~change),
         y       = expression(-log[10]~p[adj]),
         colour  = NULL,
         caption = sprintf("Up: %d  |  Down: %d", n_up, n_down)) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(size = 10))

  out_pdf <- file.path(DEG_DIR, paste0("volcano_", tag, ".pdf"))
  ggsave(out_pdf, gg, width = 5, height = 4.5)

  return(data.frame(
    cell_type  = ct,
    contrast   = paste0(cond_num, "_vs_", cond_den),
    n_up       = n_up,
    n_down     = n_down,
    n_total_sig = n_up + n_down
  ))
}


# ── Run all cell types × contrasts ────────────────────────────────────────────
cell_types <- unique(meta$cell_type)
cat(sprintf("Running DESeq2: %d cell types × %d contrasts\n\n",
            length(cell_types), length(CONTRASTS)))

summary_rows <- list()

for (ct in sort(cell_types)) {
  cat(sprintf("\n=== %s ===\n", ct))
  for (contrast in CONTRASTS) {
    row <- run_deseq2(ct, contrast[1], contrast[2])
    if (!is.null(row)) summary_rows[[length(summary_rows) + 1]] <- row
  }
}

# ── Summary table ──────────────────────────────────────────────────────────────
if (length(summary_rows) > 0) {
  summary_df <- do.call(rbind, summary_rows)
  out_summary <- file.path(DEG_DIR, "summary.csv")
  write.csv(summary_df, out_summary, row.names = FALSE)
  cat("\n===== DEG summary =====\n")
  print(summary_df, row.names = FALSE)
  cat(sprintf("\nSummary saved → %s\n", out_summary))
} else {
  cat("\nNo contrasts had sufficient replicates.\n")
}

cat("\nDone. Results in:", DEG_DIR, "\n")
