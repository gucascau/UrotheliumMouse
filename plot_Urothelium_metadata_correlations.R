#!/usr/bin/env Rscript

################################################################################
# Urothelium metadata correlation plots
#
# This script loads the saved Urothelium Harmony Seurat object and generates:
#   1. Cramer's V association heatmap between categorical metadata columns.
#   2. Cluster-composition correlation heatmaps for each metadata column.
#
# It is intentionally separate from cluster_Urothelium_harmony.R so these plots
# can be regenerated without rerunning Harmony.
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

DATA_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/"
OUT_DIR <- file.path(DATA_DIR, "integration_output")
UroDir <- file.path(OUT_DIR, "UrotheliumAll")

args <- commandArgs(trailingOnly = TRUE)
RDS_PATH <- if (length(args) >= 1) {
  args[[1]]
} else {
  file.path(OUT_DIR, "Urothelium_harmony_integrated.rds")
}

COR_DIR <- if (length(args) >= 2) {
  args[[2]]
} else {
  file.path(UroDir, "metadata_correlations")
}

METADATA_COR_COLS <- c("condition", "sample_id", "technology", "paper", "tissue")
CLUSTER_COL <- "seurat_clusters"
MIN_CELLS_PER_LEVEL <- 20
DROP_CONDITIONS <- c("E9To13.5Gestation", "E18_5_Kidney")

dir.create(COR_DIR, showWarnings = FALSE, recursive = TRUE)

safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("_+", "_", x)
}

sorted_levels <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x)]
  if (length(x) == 0) return(character(0))
  if (all(grepl("^[0-9]+$", x))) {
    return(as.character(sort(as.integer(x))))
  }
  sort(x)
}

cramers_v <- function(x, y) {
  keep <- !is.na(x) & !is.na(y)
  x <- factor(as.character(x[keep]))
  y <- factor(as.character(y[keep]))

  if (length(x) < 2 || nlevels(x) < 2 || nlevels(y) < 2) return(NA_real_)

  tab <- table(x, y)
  if (any(dim(tab) < 2)) return(NA_real_)

  chi2 <- suppressWarnings(chisq.test(tab, correct = FALSE)$statistic)
  n <- sum(tab)
  r <- nrow(tab)
  k <- ncol(tab)

  # Bias-corrected Cramer's V, useful when category counts are unbalanced.
  phi2 <- as.numeric(chi2) / n
  phi2_corr <- max(0, phi2 - ((k - 1) * (r - 1)) / (n - 1))
  r_corr <- r - ((r - 1)^2) / (n - 1)
  k_corr <- k - ((k - 1)^2) / (n - 1)
  denom <- min(k_corr - 1, r_corr - 1)

  if (is.na(denom) || denom <= 0) return(NA_real_)
  sqrt(phi2_corr / denom)
}

plot_association_heatmap <- function(mat, title, out_pdf, out_csv) {
  write.csv(mat, out_csv, quote = FALSE)

  df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(df) <- c("metadata_x", "metadata_y", "association")
  df$metadata_x <- factor(df$metadata_x, levels = colnames(mat))
  df$metadata_y <- factor(df$metadata_y, levels = rev(rownames(mat)))

  p <- ggplot(df, aes(x = metadata_x, y = metadata_y, fill = association)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(is.na(association), "", sprintf("%.2f", association))),
              size = 3.6) +
    scale_fill_gradient(
      low = "#f7fbff", high = "#084081", limits = c(0, 1),
      na.value = "grey90"
    ) +
    coord_fixed() +
    labs(title = title, x = NULL, y = NULL, fill = "Cramer's V") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )

  ggsave(out_pdf, plot = p, width = 7, height = 6)
  invisible(p)
}

plot_correlation_heatmap <- function(mat, title, out_pdf, out_csv) {
  write.csv(mat, out_csv, quote = FALSE)

  df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(df) <- c("level_x", "level_y", "correlation")
  df$level_x <- factor(df$level_x, levels = colnames(mat))
  df$level_y <- factor(df$level_y, levels = rev(rownames(mat)))

  n_levels <- ncol(mat)
  axis_text_size <- ifelse(n_levels > 45, 4.5, ifelse(n_levels > 25, 6, 8))
  show_values <- n_levels <= 20

  p <- ggplot(df, aes(x = level_x, y = level_y, fill = correlation)) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradient2(
      low = "#8c510a", mid = "white", high = "#01665e",
      midpoint = 0, limits = c(-1, 1), na.value = "grey90"
    ) +
    coord_fixed() +
    labs(title = title, x = NULL, y = NULL, fill = "Pearson r") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = axis_text_size),
      axis.text.y = element_text(size = axis_text_size),
      panel.grid = element_blank()
    )

  if (show_values) {
    p <- p +
      geom_text(aes(label = ifelse(is.na(correlation), "", sprintf("%.2f", correlation))),
                size = 3)
  }

  width <- max(6, min(32, 4 + n_levels * 0.35))
  height <- max(5, min(32, 3.5 + n_levels * 0.35))
  ggsave(out_pdf, plot = p, width = width, height = height)
  invisible(p)
}

message("Loading saved Urothelium object: ", RDS_PATH)
if (!file.exists(RDS_PATH)) {
  stop("Input RDS not found: ", RDS_PATH)
}

so <- readRDS(RDS_PATH)
meta <- so@meta.data

message("Loaded metadata for ", format(nrow(meta), big.mark = ","), " cells.")

if ("condition" %in% colnames(meta)) {
  keep_cells <- !(as.character(meta$condition) %in% DROP_CONDITIONS)
  so <- so[, keep_cells]
  meta <- so@meta.data
  message(
    "After dropping developmental conditions (",
    paste(DROP_CONDITIONS, collapse = ", "),
    "): ",
    format(nrow(meta), big.mark = ","),
    " cells."
  )
}

message("Computing average expression per group ...")
avg_expr_cols <- intersect(METADATA_COR_COLS, colnames(meta))
avg_expr_list <- setNames(lapply(avg_expr_cols, function(grp) {
  AverageExpression(so, group.by = grp, assay = "RNA")[["RNA"]]
}), avg_expr_cols)

rm(so)
gc()

metadata_cols <- intersect(METADATA_COR_COLS, colnames(meta))
if (length(metadata_cols) == 0) {
  stop("None of the requested metadata columns were found: ",
       paste(METADATA_COR_COLS, collapse = ", "))
}

metadata_summary <- data.frame(
  metadata = metadata_cols,
  n_levels = vapply(metadata_cols, function(x) {
    length(unique(as.character(meta[[x]][!is.na(meta[[x]])])))
  }, integer(1)),
  n_nonmissing_cells = vapply(metadata_cols, function(x) {
    sum(!is.na(meta[[x]]))
  }, integer(1)),
  stringsAsFactors = FALSE
)
write.csv(metadata_summary,
          file.path(COR_DIR, "metadata_columns_summary.csv"),
          row.names = FALSE, quote = FALSE)

if (length(metadata_cols) >= 2) {
  assoc_mat <- matrix(
    NA_real_,
    nrow = length(metadata_cols),
    ncol = length(metadata_cols),
    dimnames = list(metadata_cols, metadata_cols)
  )

  for (i in seq_along(metadata_cols)) {
    for (j in seq_along(metadata_cols)) {
      if (i == j) {
        assoc_mat[i, j] <- 1
      } else {
        assoc_mat[i, j] <- cramers_v(meta[[metadata_cols[i]]], meta[[metadata_cols[j]]])
      }
    }
  }

  plot_association_heatmap(
    assoc_mat,
    "Metadata association after removing developmental datasets",
    file.path(COR_DIR, "metadata_CramersV_association_heatmap.pdf"),
    file.path(COR_DIR, "metadata_CramersV_association_matrix.csv")
  )
} else {
  warning("Only one metadata column found, so Cramer's V heatmap was skipped.")
}

if (!(CLUSTER_COL %in% colnames(meta))) {
  stop("Cluster column not found: ", CLUSTER_COL)
}

for (grp in metadata_cols) {
  meta_df <- meta[, c(CLUSTER_COL, grp), drop = FALSE]
  colnames(meta_df) <- c("cluster", "group")
  meta_df <- meta_df[!is.na(meta_df$cluster) & !is.na(meta_df$group), , drop = FALSE]
  meta_df$cluster <- as.character(meta_df$cluster)
  meta_df$group <- as.character(meta_df$group)

  group_counts <- sort(table(meta_df$group), decreasing = TRUE)
  keep_groups <- names(group_counts[group_counts >= MIN_CELLS_PER_LEVEL])

  if (length(keep_groups) < 2) {
    warning("Skipping ", grp, ": fewer than two levels with at least ",
            MIN_CELLS_PER_LEVEL, " cells.")
    next
  }

  meta_df <- meta_df[meta_df$group %in% keep_groups, , drop = FALSE]
  cluster_levels <- sorted_levels(meta_df$cluster)
  group_levels <- keep_groups

  counts_mat <- table(
    cluster = factor(meta_df$cluster, levels = cluster_levels),
    group = factor(meta_df$group, levels = group_levels)
  )
  prop_mat <- prop.table(counts_mat, margin = 2)

  if (nrow(prop_mat) < 2) {
    warning("Skipping ", grp, ": cluster-composition correlation needs at least two clusters.")
    next
  }

  cor_mat <- cor(as.matrix(prop_mat), method = "pearson", use = "pairwise.complete.obs")

  prefix <- paste0("metadata_", safe_filename(grp), "_cluster_composition")
  write.csv(as.data.frame.matrix(counts_mat),
            file.path(COR_DIR, paste0(prefix, "_counts.csv")),
            quote = FALSE)
  write.csv(as.data.frame.matrix(prop_mat),
            file.path(COR_DIR, paste0(prefix, "_proportions.csv")),
            quote = FALSE)

  plot_correlation_heatmap(
    cor_mat,
    paste0(grp, ": correlation of cluster-composition profiles"),
    file.path(COR_DIR, paste0(prefix, "_correlation_heatmap.pdf")),
    file.path(COR_DIR, paste0(prefix, "_correlation_matrix.csv"))
  )
}

# ── Average expression correlation ────────────────────────────────────────────
message("Plotting average expression correlation heatmaps ...")

for (grp in metadata_cols) {
  if (!(grp %in% names(avg_expr_list))) next
  avg_mat <- avg_expr_list[[grp]]

  group_counts <- table(meta[[grp]][!is.na(meta[[grp]])])
  keep_groups <- names(group_counts[group_counts >= MIN_CELLS_PER_LEVEL])
  keep_groups <- intersect(keep_groups, colnames(avg_mat))

  if (length(keep_groups) < 2) {
    warning("Skipping avg-expr correlation for ", grp,
            ": fewer than two levels with at least ", MIN_CELLS_PER_LEVEL, " cells.")
    next
  }

  avg_mat <- avg_mat[, keep_groups, drop = FALSE]
  cor_mat <- cor(as.matrix(avg_mat), method = "pearson", use = "pairwise.complete.obs")

  prefix <- paste0("avgexpr_", safe_filename(grp))
  plot_correlation_heatmap(
    cor_mat,
    paste0(grp, ": correlation of average expression profiles"),
    file.path(COR_DIR, paste0(prefix, "_correlation_heatmap.pdf")),
    file.path(COR_DIR, paste0(prefix, "_correlation_matrix.csv"))
  )
}

generated_files <- sort(list.files(COR_DIR, full.names = FALSE))
writeLines(generated_files, file.path(COR_DIR, "metadata_correlation_files.txt"))

message("Metadata correlation outputs saved to: ", COR_DIR)
message("Generated ", length(generated_files), " files.")
