#!/usr/bin/env Rscript
# 08b_pseudobulk_deg_rds.R — Pseudobulk DEG from the Seurat RDS file.
#
# NOTE: Requires 07b_convert_annotated_to_rds.R to have been run with the
#       corrected conversion (counts="counts", data="X") so that the counts
#       slot holds raw integer counts.  The old RDS (counts="X") will give
#       wrong results with DESeq2 — re-run 07b first if needed.
#
# Limitation: only the 3,000 HVGs are tested.  For full-transcriptome DEG
#             use 08_pseudobulk.py + 09_deseq2_deg.R instead.
#
# Outputs (in output/pseudobulk_rds/DEG/)
# ----------------------------------------
#   <cell_type>_<contrast>.csv      — full DESeq2 results
#   volcano_<cell_type>_<contrast>.pdf
#   summary.csv                     — sig DEG counts per cell type × contrast
#
# Usage
# -----
#   Rscript 08b_pseudobulk_deg_rds.R

suppressPackageStartupMessages({
  library(Seurat)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})

# ── Parameters ─────────────────────────────────────────────────────────────────
BASE_DIR  <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
RDS_PATH  <- file.path(BASE_DIR, "RenalUrotheliumScripts", "output",
                       "RenalUrothelium_integrated_annotated.rds")
OUT_DIR   <- file.path(BASE_DIR, "RenalUrotheliumScripts", "output",
                       "pseudobulk_rds", "DEG")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Cell type column to use for grouping
CELL_TYPE_KEY <- "scanvi_Lake_label"   # or: scanvi_MKA_label, leiden_scVI

# Experimental sources only — exclude reference atlases
EXPERIMENTAL_SOURCES <- c(
  "KidneyHealthy1", "KidneyHealthy2", "KidneyHealthy3",
  "KidneyHealthy4", "KidneyHealthy5",
  "KidneyTET2UUO",
  "KidneyUUO1", "KidneyUUO2", "KidneyUUO3", "KidneyUUO4",
  "KidneyUUO5", "KidneyUUO6", "KidneyUUO7", "KidneyUUO8",
  "KidneyrUUO1"
)

# Contrasts: each row is c(numerator, denominator)
CONTRASTS <- list(
  c("UUO",     "Healthy"),
  c("rUUO",    "Healthy"),
  c("TET2UUO", "Healthy"),
  c("UUO",     "rUUO")
)

# Minimum cells per pseudobulk group; minimum replicates per condition
MIN_CELLS      <- 10
MIN_REPLICATES <- 3

# DEG thresholds
PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1.0    # |log2FC|

# ── Load and subset Seurat object ──────────────────────────────────────────────
cat("Loading Seurat RDS ...\n")
so <- readRDS(RDS_PATH)
cat(sprintf("  %d cells × %d genes\n", ncol(so), nrow(so)))

# Verify counts slot has integer-like values (sanity check)
counts_sample <- GetAssayData(so, slot = "counts")[1:5, 1:5]
if (any(counts_sample != round(counts_sample))) {
  stop(
    "counts slot contains non-integer values — the RDS was built with the old\n",
    "07b script (counts='X'). Please re-run 07b_convert_annotated_to_rds.R\n",
    "with counts='counts', data='X', then rerun this script."
  )
}

# Subset to experimental samples
so <- so[, so@meta.data$source %in% EXPERIMENTAL_SOURCES]
cat(sprintf("  After experimental filter: %d cells\n", ncol(so)))
cat("  Conditions:", paste(unique(so@meta.data$condition), collapse = ", "), "\n\n")

# ── Pseudobulk aggregation ─────────────────────────────────────────────────────
# Sum raw counts per sample_id × cell_type group
build_pseudobulk <- function(so, cell_type_col) {
  meta    <- so@meta.data
  counts  <- GetAssayData(so, slot = "counts")   # genes × cells

  groups  <- paste0(meta$sample_id, "__", meta[[cell_type_col]])
  unique_groups <- unique(groups)

  pb_list <- lapply(unique_groups, function(g) {
    idx <- which(groups == g)
    if (length(idx) < MIN_CELLS) return(NULL)
    Matrix::rowSums(counts[, idx, drop = FALSE])
  })
  names(pb_list) <- unique_groups

  # Drop groups below MIN_CELLS
  pb_list <- Filter(Negate(is.null), pb_list)

  pb_counts <- do.call(cbind, pb_list)
  colnames(pb_counts) <- names(pb_list)

  # Build metadata for each pseudobulk sample
  pb_meta <- do.call(rbind, lapply(names(pb_list), function(g) {
    idx <- which(groups == g)[1]
    data.frame(
      pseudobulk_id = g,
      sample_id     = meta$sample_id[idx],
      condition     = meta$condition[idx],
      cell_type     = meta[[cell_type_col]][idx],
      n_cells       = sum(groups == g),
      row.names     = g,
      stringsAsFactors = FALSE
    )
  }))

  list(counts = pb_counts, meta = pb_meta)
}

cat("Building pseudobulk count matrix ...\n")
pb      <- build_pseudobulk(so, CELL_TYPE_KEY)
pb_counts <- pb$counts
pb_meta   <- pb$meta

cat(sprintf("  %d genes × %d pseudobulk samples\n",
            nrow(pb_counts), ncol(pb_counts)))
cat("  Cell types:\n")
print(table(pb_meta$cell_type))

# ── DESeq2 per cell type × contrast ───────────────────────────────────────────
run_deseq2 <- function(ct, cond_num, cond_den) {
  tag <- paste0(gsub("[/ ]", "_", ct), "__", cond_num, "_vs_", cond_den)

  sel      <- pb_meta$cell_type == ct & pb_meta$condition %in% c(cond_num, cond_den)
  meta_sub <- pb_meta[sel, , drop = FALSE]
  cnt_sub  <- pb_counts[, rownames(meta_sub), drop = FALSE]

  n_num <- sum(meta_sub$condition == cond_num)
  n_den <- sum(meta_sub$condition == cond_den)

  if (n_num < MIN_REPLICATES || n_den < MIN_REPLICATES) {
    cat(sprintf("  [%s] SKIP — %s n=%d, %s n=%d (need ≥%d)\n",
                tag, cond_num, n_num, cond_den, n_den, MIN_REPLICATES))
    return(NULL)
  }

  cat(sprintf("  [%s] %s(n=%d) vs %s(n=%d) ...\n",
              ct, cond_num, n_num, cond_den, n_den))

  # Low-count filter: ≥10 counts in ≥2 pseudobulk samples
  keep    <- rowSums(cnt_sub >= 10) >= 2
  cnt_sub <- cnt_sub[keep, , drop = FALSE]
  cat(sprintf("    Genes after filter: %d / %d\n", sum(keep), length(keep)))

  meta_sub$condition <- factor(meta_sub$condition,
                               levels = c(cond_den, cond_num))

  dds <- tryCatch(
    DESeqDataSetFromMatrix(
      countData = as.matrix(cnt_sub),
      colData   = meta_sub,
      design    = ~ condition
    ),
    error = function(e) { cat("    DESeqDataSetFromMatrix:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(dds)) return(NULL)

  dds <- tryCatch(
    DESeq(dds, quiet = TRUE),
    error = function(e) { cat("    DESeq():", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(dds)) return(NULL)

  res <- results(dds,
                 contrast = c("condition", cond_num, cond_den),
                 alpha    = PADJ_CUTOFF)
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    arrange(padj, desc(abs(log2FoldChange)))

  write.csv(res_df, file.path(OUT_DIR, paste0(tag, ".csv")), row.names = FALSE)

  n_up   <- sum(res_df$padj < PADJ_CUTOFF & res_df$log2FoldChange >  LFC_CUTOFF, na.rm = TRUE)
  n_down <- sum(res_df$padj < PADJ_CUTOFF & res_df$log2FoldChange < -LFC_CUTOFF, na.rm = TRUE)
  cat(sprintf("    Up: %d  Down: %d\n", n_up, n_down))

  # Volcano plot
  plot_df <- res_df %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(
      direction = case_when(
        padj < PADJ_CUTOFF & log2FoldChange >  LFC_CUTOFF ~ "Up",
        padj < PADJ_CUTOFF & log2FoldChange < -LFC_CUTOFF ~ "Down",
        TRUE ~ "NS"
      ),
      label = ifelse(
        padj < PADJ_CUTOFF & abs(log2FoldChange) >= 2 &
          rank(padj) <= 20, gene, NA
      )
    )

  gg <- ggplot(plot_df, aes(log2FoldChange, -log10(padj), colour = direction)) +
    geom_point(size = 0.7, alpha = 0.6) +
    scale_colour_manual(values = c(Up = "#d62728", Down = "#1f77b4", NS = "grey70")) +
    geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF),
               linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    geom_hline(yintercept = -log10(PADJ_CUTOFF),
               linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    geom_text_repel(aes(label = label), size = 2.5,
                    max.overlaps = 20, na.rm = TRUE) +
    labs(
      title   = paste0(ct, "  |  ", cond_num, " vs ", cond_den,
                       "  [", CELL_TYPE_KEY, "]"),
      x       = expression(log[2]~fold~change),
      y       = expression(-log[10]~p[adj]),
      colour  = NULL,
      caption = paste0("Up: ", n_up, "  |  Down: ", n_down,
                       "  (padj<", PADJ_CUTOFF, ", |LFC|>", LFC_CUTOFF, ")\n",
                       "Note: 3,000 HVGs only — use 08_pseudobulk.py for full transcriptome")
    ) +
    theme_classic(base_size = 11)

  ggsave(file.path(OUT_DIR, paste0("volcano_", tag, ".pdf")),
         gg, width = 5.5, height = 5)

  data.frame(cell_type   = ct,
             contrast    = paste0(cond_num, "_vs_", cond_den),
             n_genes_tested = sum(keep),
             n_up        = n_up,
             n_down      = n_down)
}

cell_types   <- sort(unique(pb_meta$cell_type))
summary_rows <- list()

cat(sprintf("\nRunning DESeq2: %d cell types × %d contrasts ...\n",
            length(cell_types), length(CONTRASTS)))

for (ct in cell_types) {
  cat(sprintf("\n=== %s ===\n", ct))
  for (contrast in CONTRASTS) {
    row <- run_deseq2(ct, contrast[1], contrast[2])
    if (!is.null(row)) summary_rows[[length(summary_rows) + 1]] <- row
  }
}

if (length(summary_rows) > 0) {
  summary_df <- do.call(rbind, summary_rows)
  write.csv(summary_df, file.path(OUT_DIR, "summary.csv"), row.names = FALSE)
  cat("\n===== DEG summary =====\n")
  print(summary_df, row.names = FALSE)
}

cat("\nDone. Results in:", OUT_DIR, "\n")
cat("Reminder: results cover only the 3,000 HVGs.\n",
    "For full-transcriptome DEG run 08_pseudobulk.py + 09_deseq2_deg.R\n")
