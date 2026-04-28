#!/usr/bin/env Rscript

################################################################################
# 06_metadata_correlations.R
#
# Input : output/RenalUrothelium_uro_cells_fullgene_harmony_integrated.rds
#
# Outputs (in output/metadata_correlations/):
#   1. Cramer's V association heatmap between categorical metadata columns
#   2. Average-expression Pearson correlation heatmaps per metadata column
#      (genes x groups average log-norm expression → group-group Pearson r)
#   3. Supporting CSVs (average expression matrices, correlation matrices)
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(Matrix)
})

BASE_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
OUT_DIR  <- file.path(BASE_DIR, "RenalUrotheliumScripts", "output")
COR_DIR  <- file.path(OUT_DIR, "metadata_correlations")

args     <- commandArgs(trailingOnly = TRUE)
RDS_PATH <- if (length(args) >= 1) args[[1]] else
  file.path(OUT_DIR, "RenalUrothelium_uro_cells_fullgene_harmony_integrated.rds")
if (length(args) >= 2) COR_DIR <- args[[2]]

dir.create(COR_DIR, showWarnings = FALSE, recursive = TRUE)

METADATA_COR_COLS   <- c("condition", "sample_id", "technology", "paper", "tissue")
MIN_CELLS_PER_LEVEL <- 20
MIN_EXPR_MEAN       <- 0.01   # filter out genes with very low average expression
TOP_VAR_GENES       <- 2000   # use top variable genes for correlation

# ── Helpers ───────────────────────────────────────────────────────────────────

safe_filename <- function(x) {
  gsub("_+", "_", gsub("[^A-Za-z0-9]+", "_", x))
}

cramers_v <- function(x, y) {
  keep <- !is.na(x) & !is.na(y)
  x <- factor(as.character(x[keep]))
  y <- factor(as.character(y[keep]))
  if (length(x) < 2 || nlevels(x) < 2 || nlevels(y) < 2) return(NA_real_)
  tab <- table(x, y)
  if (any(dim(tab) < 2)) return(NA_real_)
  chi2  <- suppressWarnings(chisq.test(tab, correct = FALSE)$statistic)
  n     <- sum(tab)
  r     <- nrow(tab)
  k     <- ncol(tab)
  phi2c <- max(0, as.numeric(chi2) / n - ((k - 1) * (r - 1)) / (n - 1))
  rc    <- r - ((r - 1)^2) / (n - 1)
  kc    <- k - ((k - 1)^2) / (n - 1)
  denom <- min(kc - 1, rc - 1)
  if (is.na(denom) || denom <= 0) return(NA_real_)
  sqrt(phi2c / denom)
}

plot_association_heatmap <- function(mat, title, out_pdf, out_csv) {
  write.csv(mat, out_csv, quote = FALSE)
  df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(df) <- c("metadata_x", "metadata_y", "association")
  df$metadata_x <- factor(df$metadata_x, levels = colnames(mat))
  df$metadata_y <- factor(df$metadata_y, levels = rev(rownames(mat)))
  p <- ggplot(df, aes(x = metadata_x, y = metadata_y, fill = association)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(is.na(association), "",
                                 sprintf("%.2f", association))), size = 3.6) +
    scale_fill_gradient(low = "#f7fbff", high = "#084081",
                        limits = c(0, 1), na.value = "grey90") +
    coord_fixed() +
    labs(title = title, x = NULL, y = NULL, fill = "Cramer's V") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank())
  ggsave(out_pdf, plot = p, width = 7, height = 6)
  invisible(p)
}

plot_correlation_heatmap <- function(mat, title, out_pdf, out_csv) {
  write.csv(mat, out_csv, quote = FALSE)
  df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(df) <- c("level_x", "level_y", "correlation")
  df$level_x <- factor(df$level_x, levels = colnames(mat))
  df$level_y <- factor(df$level_y, levels = rev(rownames(mat)))
  n   <- ncol(mat)
  txt <- ifelse(n > 45, 4.5, ifelse(n > 25, 6, 8))
  p <- ggplot(df, aes(x = level_x, y = level_y, fill = correlation)) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradient2(low = "#8c510a", mid = "white", high = "#01665e",
                         midpoint = 0, limits = c(-1, 1), na.value = "grey90") +
    coord_fixed() +
    labs(title = title, x = NULL, y = NULL, fill = "Pearson r") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = txt),
          axis.text.y = element_text(size = txt),
          panel.grid = element_blank())
  if (n <= 30)
    p <- p + geom_text(aes(label = ifelse(is.na(correlation), "",
                                          sprintf("%.2f", correlation))), size = 3)
  w <- max(6, min(32, 4 + n * 0.35))
  h <- max(5, min(32, 3.5 + n * 0.35))
  ggsave(out_pdf, plot = p, width = w, height = h)
  invisible(p)
}

# ── Load ──────────────────────────────────────────────────────────────────────

message("Loading: ", RDS_PATH)
if (!file.exists(RDS_PATH)) stop("Input RDS not found: ", RDS_PATH)

so <- readRDS(RDS_PATH)
message(sprintf("  Loaded: %d cells × %d genes", ncol(so), nrow(so)))

# ── Column inventory ──────────────────────────────────────────────────────────

meta          <- so@meta.data
metadata_cols <- intersect(METADATA_COR_COLS, colnames(meta))
if (length(metadata_cols) == 0)
  stop("None of the requested metadata columns found: ",
       paste(METADATA_COR_COLS, collapse = ", "))
message("  Metadata columns: ", paste(metadata_cols, collapse = ", "))

summary_df <- data.frame(
  metadata     = metadata_cols,
  n_levels     = vapply(metadata_cols, function(x)
    length(unique(as.character(meta[[x]][!is.na(meta[[x]])]))), integer(1)),
  n_nonmissing = vapply(metadata_cols, function(x)
    sum(!is.na(meta[[x]])), integer(1)),
  stringsAsFactors = FALSE
)
write.csv(summary_df, file.path(COR_DIR, "metadata_columns_summary.csv"),
          row.names = FALSE, quote = FALSE)
print(summary_df)

# ── 1. Cramer's V association heatmap ────────────────────────────────────────

if (length(metadata_cols) >= 2) {
  message("Computing Cramer's V association matrix ...")
  assoc_mat <- matrix(NA_real_, nrow = length(metadata_cols),
                      ncol = length(metadata_cols),
                      dimnames = list(metadata_cols, metadata_cols))
  for (i in seq_along(metadata_cols))
    for (j in seq_along(metadata_cols))
      assoc_mat[i, j] <- if (i == j) 1 else
        cramers_v(meta[[metadata_cols[i]]], meta[[metadata_cols[j]]])

  plot_association_heatmap(
    assoc_mat,
    "Metadata association (renal urothelial cells)",
    file.path(COR_DIR, "metadata_CramersV_association_heatmap.pdf"),
    file.path(COR_DIR, "metadata_CramersV_association_matrix.csv")
  )
  message("  Saved: metadata_CramersV_association_heatmap.pdf")
} else {
  warning("Only one metadata column — Cramer's V heatmap skipped.")
}

# ── 2. Average-expression Pearson correlation heatmaps ───────────────────────
#
# Strategy:
#   a. Pull log-normalised data layer (genes × cells)
#   b. Filter to top variable genes by row variance across the full matrix
#   c. For each grouping variable, compute per-group column means
#      → genes × groups average expression matrix
#   d. Correlate groups (columns) by Pearson r → groups × groups matrix
#   e. Plot + save

message(sprintf(
  "Extracting log-norm data layer (top %d variable genes, mean >= %.3f) ...",
  TOP_VAR_GENES, MIN_EXPR_MEAN))

expr_mat <- GetAssayData(so, assay = "RNA", layer = "data")  # genes × cells (sparse)

# Filter to genes expressed above minimum mean
gene_means <- Matrix::rowMeans(expr_mat)
expr_mat   <- expr_mat[gene_means >= MIN_EXPR_MEAN, ]
message(sprintf("  Genes after mean filter: %d", nrow(expr_mat)))

# Select top variable genes by row variance (sparse-safe)
gene_vars  <- sparseMatrixStats::rowVars(expr_mat)
if (!requireNamespace("sparseMatrixStats", quietly = TRUE)) {
  # fallback: row variance via rowMeans on dense chunk if sparseMatrixStats absent
  gene_means2 <- Matrix::rowMeans(expr_mat)
  gene_vars   <- Matrix::rowMeans(expr_mat^2) - gene_means2^2
}
top_genes  <- order(gene_vars, decreasing = TRUE)[seq_len(min(TOP_VAR_GENES, length(gene_vars)))]
expr_mat   <- expr_mat[top_genes, ]
message(sprintf("  Using %d variable genes for correlation.", nrow(expr_mat)))

for (grp in metadata_cols) {
  grp_vec <- as.character(meta[[grp]])
  grp_vec[is.na(grp_vec)] <- NA_character_

  grp_counts  <- sort(table(grp_vec), decreasing = TRUE)
  keep_groups <- names(grp_counts[grp_counts >= MIN_CELLS_PER_LEVEL])

  if (length(keep_groups) < 2) {
    warning(sprintf("Skipping %s: fewer than 2 levels with >= %d cells.",
                    grp, MIN_CELLS_PER_LEVEL))
    next
  }

  keep_cells <- grp_vec %in% keep_groups
  sub_expr   <- expr_mat[, keep_cells, drop = FALSE]
  sub_grp    <- grp_vec[keep_cells]

  # Average expression per group: genes × groups
  avg_mat <- do.call(cbind, lapply(keep_groups, function(g) {
    idx <- sub_grp == g
    Matrix::rowMeans(sub_expr[, idx, drop = FALSE])
  }))
  colnames(avg_mat) <- keep_groups
  rownames(avg_mat) <- rownames(sub_expr)

  # Pearson correlation between groups based on their avg expression profiles
  cor_mat <- cor(avg_mat, method = "pearson", use = "pairwise.complete.obs")

  prefix <- paste0("avgexpr_", safe_filename(grp))
  write.csv(as.data.frame(avg_mat),
            file.path(COR_DIR, paste0(prefix, "_average_expression.csv")),
            quote = FALSE)
  plot_correlation_heatmap(
    cor_mat,
    paste0(grp, ": avg-expression correlation (renal urothelial cells)"),
    file.path(COR_DIR, paste0(prefix, "_correlation_heatmap.pdf")),
    file.path(COR_DIR, paste0(prefix, "_correlation_matrix.csv"))
  )
  message(sprintf("  Saved: %s_correlation_heatmap.pdf", prefix))
}

# ── Summary ───────────────────────────────────────────────────────────────────

generated <- sort(list.files(COR_DIR, full.names = FALSE))
writeLines(generated, file.path(COR_DIR, "metadata_correlation_files.txt"))
message(sprintf("\nDone — %d files saved to: %s", length(generated), COR_DIR))
