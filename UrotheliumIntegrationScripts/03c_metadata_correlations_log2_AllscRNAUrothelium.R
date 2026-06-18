################################################################################
# 03c_metadata_correlations_log2_AllscRNAUrothelium.R
#
# Improved version of 03b: instead of correlating raw rowMeans of log1p-
# normalised data, this script:
#   1. Back-transforms log1p → linear normalised counts (expm1)
#   2. Computes per-group arithmetic mean on the linear scale
#   3. Applies log2(avg + 1) so the final matrix is on a log2 scale
#   4. Filters near-zero genes on the log2 scale (stricter than before)
#   5. Computes BOTH Pearson and Spearman correlations for every heatmap
#
# Input : output/AllUrothelium_gated_UMAP.rds   (or first CLI arg)
# Output: output/AllUrothelium_metadata_correlations_log2/   (or second CLI arg)
#   - Cramer's V association heatmap (unchanged)
#   - Cluster-composition Pearson + Spearman heatmaps
#   - Average-expression Pearson + Spearman heatmaps (per gene mode)
#   - Supporting CSVs
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(Matrix)
  library(dplyr)
  library(tidyr)
})

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/FinalUrotheliumCells"
SCR_DIR  <- file.path(BASE_DIR, "UrotheliumScripts")
OUT_DIR  <- file.path(SCR_DIR, "output")
COR_DIR  <- file.path(OUT_DIR, "AllUrothelium_metadata_correlations_log2")

args     <- commandArgs(trailingOnly = TRUE)
RDS_PATH <- if (length(args) >= 1) args[[1]] else
  file.path(OUT_DIR, "AllUrothelium_gated_UMAP.rds")
if (length(args) >= 2) COR_DIR <- args[[2]]

dir.create(COR_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Parameters ─────────────────────────────────────────────────────────────────
METADATA_COR_COLS   <- c("LabelClass","Finalgsm_id","FinalConditionL1",
                          "FinalConditionL2","sample_id","technology",
                          "paper","tissue","Categories","Finalpaper")
CLUSTER_COL         <- "seurat_clusters"
MIN_CELLS_PER_LEVEL <- 20

# Log2-scale filtering (applied after log2(avg + 1) transformation)
LOG2_PSEUDOCOUNT      <- 1     # added before log2 to avoid log(0)
MIN_LOG2_EXPR         <- 0.5   # gene log2-mean must be >= 0.5 (≈ norm-count 0.41)
MIN_EXPR_FRAC_SAMPLES <- 0.25  # gene must be above threshold in >= 25% of groups

TOP_VAR_GENES <- 5000
DROP_CONDITIONS <- c("E9To13.5Gestation", "E18_5_Kidney")

# Gene-selection modes (same as 03b)
GENE_MODES <- c("hvg", "allExpr_noMtRibo", "allExpr_noMt")
GENE_MODE_LABELS <- c(
  hvg              = "HVG (no mt/ribo)",
  allExpr_noMtRibo = "All-expr (no mt, no ribo)",
  allExpr_noMt     = "All-expr (no mt)"
)

# ── Helpers ────────────────────────────────────────────────────────────────────

safe_filename <- function(x) {
  gsub("_+", "_", gsub("[^A-Za-z0-9]+", "_", x))
}

sorted_levels <- function(x) {
  x <- unique(as.character(x[!is.na(x)]))
  if (length(x) == 0) return(character(0))
  if (all(grepl("^[0-9]+$", x))) return(as.character(sort(as.integer(x))))
  sort(x)
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
          panel.grid  = element_blank())
  ggsave(out_pdf, plot = p, width = 7, height = 6)
  invisible(p)
}

# method_label: "Pearson r" or "Spearman ρ"
plot_correlation_heatmap <- function(mat, title, out_pdf, out_csv,
                                     method_label = "Pearson r") {
  mat_for_clust <- mat
  mat_for_clust[is.na(mat_for_clust)] <- 0
  hc  <- hclust(as.dist(1 - mat_for_clust), method = "average")
  ord <- hc$labels[hc$order]
  mat <- mat[ord, ord]

  write.csv(mat, out_csv, quote = FALSE)
  df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(df) <- c("level_x", "level_y", "correlation")
  df$level_x <- factor(df$level_x, levels = colnames(mat))
  df$level_y <- factor(df$level_y, levels = rev(rownames(mat)))

  n_levels       <- ncol(mat)
  axis_text_size <- ifelse(n_levels > 45, 4.5, ifelse(n_levels > 25, 6, 8))
  show_values    <- n_levels <= 20

  p <- ggplot(df, aes(x = level_x, y = level_y, fill = correlation)) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradient2(low = "#01665e", mid = "white", high = "#8c510a",
                         midpoint = 0, limits = c(-1, 1), na.value = "grey90") +
    coord_fixed() +
    labs(title = title, x = NULL, y = NULL, fill = method_label) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1,
                                     size = axis_text_size),
          axis.text.y = element_text(size = axis_text_size),
          panel.grid  = element_blank())

  if (show_values)
    p <- p + geom_text(aes(label = ifelse(is.na(correlation), "",
                                          sprintf("%.2f", correlation))),
                       size = 3)

  width  <- max(6, min(32, 4 + n_levels * 0.35))
  height <- max(5, min(32, 3.5 + n_levels * 0.35))
  ggsave(out_pdf, plot = p, width = width, height = height)
  invisible(p)
}

# Save both Pearson and Spearman heatmaps from the same avg_mat
save_cor_pair <- function(avg_mat, grp, mode_dir, prefix, title_base, mode_label) {
  cor_pearson  <- cor(as.matrix(avg_mat), method = "pearson",
                      use = "pairwise.complete.obs")
  cor_spearman <- cor(as.matrix(avg_mat), method = "spearman",
                      use = "pairwise.complete.obs")

  plot_correlation_heatmap(
    cor_pearson,
    sprintf("%s: avg-expression Pearson\n[%s, log2(avg+1)]", grp, mode_label),
    file.path(mode_dir, paste0(prefix, "_pearson_heatmap.pdf")),
    file.path(mode_dir, paste0(prefix, "_pearson_matrix.csv")),
    method_label = "Pearson r"
  )
  plot_correlation_heatmap(
    cor_spearman,
    sprintf("%s: avg-expression Spearman\n[%s, log2(avg+1)]", grp, mode_label),
    file.path(mode_dir, paste0(prefix, "_spearman_heatmap.pdf")),
    file.path(mode_dir, paste0(prefix, "_spearman_matrix.csv")),
    method_label = "Spearman ρ"
  )
}


################################################################################
# Load
################################################################################

message("Loading: ", RDS_PATH)
if (!file.exists(RDS_PATH))
  stop("Input RDS not found: ", RDS_PATH)

so   <- readRDS(RDS_PATH)
meta <- so@meta.data
message(sprintf("Loaded %s cells.", format(nrow(meta), big.mark = ",")))

if ("condition" %in% colnames(meta)) {
  keep <- !(as.character(meta$condition) %in% DROP_CONDITIONS)
  if (any(!keep)) {
    so   <- so[, keep]
    meta <- so@meta.data
    message(sprintf("After dropping developmental conditions: %s cells.",
                    format(nrow(meta), big.mark = ",")))
  }
}

metadata_cols <- intersect(METADATA_COR_COLS, colnames(meta))
if (length(metadata_cols) == 0)
  stop("None of the requested metadata columns found: ",
       paste(METADATA_COR_COLS, collapse = ", "))
if (!(CLUSTER_COL %in% colnames(meta)))
  stop("Cluster column not found: ", CLUSTER_COL)

metadata_summary <- data.frame(
  metadata = metadata_cols,
  n_levels = vapply(metadata_cols, function(x)
    length(unique(as.character(meta[[x]][!is.na(meta[[x]])]))), integer(1)),
  n_cells  = vapply(metadata_cols, function(x)
    sum(!is.na(meta[[x]])), integer(1)),
  stringsAsFactors = FALSE
)
write.csv(metadata_summary,
          file.path(COR_DIR, "metadata_columns_summary.csv"),
          row.names = FALSE, quote = FALSE)
message("Metadata columns used:")
print(metadata_summary)


################################################################################
# Average expression — log2(avg_normcount + 1)
#
# Strategy:
#   data layer holds log1p(normalised_count).
#   expm1() reverses log1p → normalised count on a linear scale.
#   rowMeans per group on the linear scale gives the arithmetic mean.
#   log2(avg + LOG2_PSEUDOCOUNT) puts the result on log2 scale.
#
# Why not average on log scale?
#   mean(log(x)) ≠ log(mean(x)).  Averaging on the linear scale, then
#   log-transforming, preserves the interpretable "average expression level"
#   while still compressing the dynamic range for correlation analysis.
################################################################################

message("Computing log2(avg_normcount + 1) per metadata group ...")

rna_layers <- Layers(so[["RNA"]])
if (!"data" %in% rna_layers) so <- JoinLayers(so)

all_genes <- rownames(so)
mt_genes  <- grep("^mt-",    all_genes, ignore.case = TRUE, value = TRUE)
rib_genes <- grep("^Rp[sl]", all_genes, ignore.case = TRUE, value = TRUE)

hvg_genes <- VariableFeatures(so)
if (length(hvg_genes) == 0) {
  message("  No HVGs found — using all genes for HVG mode.")
  hvg_genes <- all_genes
}
hvg_genes <- head(setdiff(hvg_genes, c(mt_genes, rib_genes)), TOP_VAR_GENES)

gene_sets <- list(
  hvg              = hvg_genes,
  allExpr_noMtRibo = setdiff(all_genes, c(mt_genes, rib_genes)),
  allExpr_noMt     = setdiff(all_genes, mt_genes)
)
message(sprintf("  Gene sets — hvg: %d | allExpr_noMtRibo: %d | allExpr_noMt: %d",
                length(gene_sets$hvg),
                length(gene_sets$allExpr_noMtRibo),
                length(gene_sets$allExpr_noMt)))

# Load log1p data once; back-transform to linear normalised counts
data_mat_full <- GetAssayData(so, assay = "RNA", layer = "data")
data_mat_full <- data_mat_full[
  intersect(gene_sets$allExpr_noMt, rownames(data_mat_full)), , drop = FALSE
]

# expm1(x) = exp(x) - 1, exactly reversing log1p; sparse-safe (expm1(0) = 0)
message("  Back-transforming log1p → linear normalised counts ...")
norm_mat_full <- expm1(data_mat_full)

compute_log2_avg_expr <- function(norm_mat, meta_obj) {
  setNames(
    lapply(metadata_cols, function(grp) {
      tryCatch({
        groups     <- as.character(meta_obj[[grp]])
        valid_grps <- names(which(table(groups[!is.na(groups)]) >= 1))
        if (length(valid_grps) < 2) return(NULL)
        avg_cols <- lapply(valid_grps, function(g) {
          idx     <- which(!is.na(groups) & groups == g)
          avg_lin <- Matrix::rowMeans(norm_mat[, idx, drop = FALSE])
          log2(avg_lin + LOG2_PSEUDOCOUNT)   # log2(avg_normcount + 1)
        })
        names(avg_cols) <- valid_grps
        do.call(cbind, avg_cols)
      }, error = function(e) {
        warning("log2-avg-expr failed for ", grp, ": ", conditionMessage(e))
        NULL
      })
    }),
    metadata_cols
  )
}

avg_expr_by_mode <- setNames(
  lapply(GENE_MODES, function(mode) {
    genes    <- intersect(gene_sets[[mode]], rownames(norm_mat_full))
    norm_mat <- norm_mat_full[genes, , drop = FALSE]
    message(sprintf("  Mode '%s': %d genes", mode, nrow(norm_mat)))
    Filter(Negate(is.null), compute_log2_avg_expr(norm_mat, so@meta.data))
  }),
  GENE_MODES
)

rm(data_mat_full, norm_mat_full, so); gc()
message("  Log2 average expression computed; Seurat object freed.")


################################################################################
# 1. Cramer's V association heatmap
################################################################################

if (length(metadata_cols) >= 2) {
  message("Computing Cramer's V association matrix ...")
  assoc_mat <- matrix(
    NA_real_,
    nrow = length(metadata_cols), ncol = length(metadata_cols),
    dimnames = list(metadata_cols, metadata_cols)
  )
  for (i in seq_along(metadata_cols))
    for (j in seq_along(metadata_cols))
      assoc_mat[i, j] <- if (i == j) 1 else
        cramers_v(meta[[metadata_cols[i]]], meta[[metadata_cols[j]]])

  plot_association_heatmap(
    assoc_mat,
    "Metadata association (Cramer's V)",
    file.path(COR_DIR, "AllUrothelium_CramersV_association_heatmap.pdf"),
    file.path(COR_DIR, "AllUrothelium_CramersV_association_matrix.csv")
  )
  message("  Saved: AllUrothelium_CramersV_association_heatmap.pdf")
} else {
  warning("Only one metadata column available — Cramer's V heatmap skipped.")
}


################################################################################
# 2. Cluster-composition correlation heatmaps (Pearson + Spearman)
################################################################################

message("Generating cluster-composition correlation heatmaps ...")

for (grp in metadata_cols) {
  meta_df <- meta[, c(CLUSTER_COL, grp), drop = FALSE]
  colnames(meta_df) <- c("cluster", "group")
  meta_df <- meta_df[!is.na(meta_df$cluster) & !is.na(meta_df$group), , drop = FALSE]
  meta_df$cluster <- as.character(meta_df$cluster)
  meta_df$group   <- as.character(meta_df$group)

  grp_counts  <- sort(table(meta_df$group), decreasing = TRUE)
  keep_groups <- names(grp_counts[grp_counts >= MIN_CELLS_PER_LEVEL])
  if (length(keep_groups) < 2) {
    warning("Skipping cluster-composition for ", grp,
            ": fewer than two levels with >= ", MIN_CELLS_PER_LEVEL, " cells.")
    next
  }

  meta_df     <- meta_df[meta_df$group %in% keep_groups, , drop = FALSE]
  clust_levels <- sorted_levels(meta_df$cluster)
  prop_mat     <- prop.table(
    table(cluster = factor(meta_df$cluster, levels = clust_levels),
          group   = factor(meta_df$group,   levels = keep_groups)),
    margin = 2
  )
  if (nrow(prop_mat) < 2) next

  prefix <- paste0("AllUrothelium_clustcomp_", safe_filename(grp))

  cor_pearson  <- cor(as.matrix(prop_mat), method = "pearson",
                      use = "pairwise.complete.obs")
  cor_spearman <- cor(as.matrix(prop_mat), method = "spearman",
                      use = "pairwise.complete.obs")

  plot_correlation_heatmap(
    cor_pearson,
    sprintf("%s: cluster-composition Pearson", grp),
    file.path(COR_DIR, paste0(prefix, "_pearson_heatmap.pdf")),
    file.path(COR_DIR, paste0(prefix, "_pearson_matrix.csv")),
    method_label = "Pearson r"
  )
  plot_correlation_heatmap(
    cor_spearman,
    sprintf("%s: cluster-composition Spearman", grp),
    file.path(COR_DIR, paste0(prefix, "_spearman_heatmap.pdf")),
    file.path(COR_DIR, paste0(prefix, "_spearman_matrix.csv")),
    method_label = "Spearman ρ"
  )
  message(sprintf("  Saved Pearson + Spearman: %s", prefix))
}


################################################################################
# 3. Average-expression correlation heatmaps — log2(avg+1), Pearson + Spearman
#
# Filtering (on log2 scale, applied before correlation):
#   a. row mean across groups >= MIN_LOG2_EXPR (removes near-zero genes)
#   b. gene expressed (>= MIN_LOG2_EXPR) in >= MIN_EXPR_FRAC_SAMPLES fraction
#      of groups (removes structural zeros / tissue-specific genes)
################################################################################

for (mode in GENE_MODES) {
  mode_dir <- file.path(COR_DIR, mode)
  dir.create(mode_dir, showWarnings = FALSE, recursive = TRUE)
  message(sprintf(
    "\nAverage-expression (log2) correlations [%s] ...",
    GENE_MODE_LABELS[[mode]]
  ))

  avg_expr_list <- avg_expr_by_mode[[mode]]

  for (grp in names(avg_expr_list)) {
    avg_mat <- avg_expr_list[[grp]]

    # Keep only groups with enough cells
    grp_counts  <- table(meta[[grp]][!is.na(meta[[grp]])])
    keep_groups <- names(grp_counts[grp_counts >= MIN_CELLS_PER_LEVEL])
    keep_groups <- intersect(keep_groups, colnames(avg_mat))
    if (length(keep_groups) < 2) {
      warning("Skipping avg-expr for ", grp, " (mode=", mode,
              "): fewer than two qualifying levels.")
      next
    }
    avg_mat <- avg_mat[, keep_groups, drop = FALSE]

    # (a) Filter near-zero genes: log2 row mean < MIN_LOG2_EXPR
    gene_log2_means <- rowMeans(as.matrix(avg_mat))
    avg_mat <- avg_mat[gene_log2_means >= MIN_LOG2_EXPR, , drop = FALSE]
    message(sprintf("    [%s | %s] after mean filter: %d genes",
                    mode, grp, nrow(avg_mat)))

    # (b) Filter genes not expressed above threshold in enough groups
    n_expressed <- rowSums(as.matrix(avg_mat) >= MIN_LOG2_EXPR)
    min_samples <- ceiling(MIN_EXPR_FRAC_SAMPLES * ncol(avg_mat))
    avg_mat <- avg_mat[n_expressed >= min_samples, , drop = FALSE]
    message(sprintf("    [%s | %s] after breadth filter (>= %d groups): %d genes",
                    mode, grp, min_samples, nrow(avg_mat)))

    if (nrow(avg_mat) < 10) {
      warning("Skipping ", grp, " (mode=", mode,
              "): fewer than 10 genes remain after filtering.")
      next
    }

    prefix <- paste0("AllUrothelium_avgexpr_log2_", safe_filename(grp))

    # Save filtered log2 matrix
    write.csv(avg_mat,
              file.path(mode_dir, paste0(prefix, "_log2_avg_expression.csv")),
              quote = FALSE)

    # Pearson and Spearman heatmaps
    save_cor_pair(avg_mat, grp, mode_dir, prefix,
                  title_base = grp, mode_label = GENE_MODE_LABELS[[mode]])
    message(sprintf("  Saved Pearson + Spearman: %s/%s", mode, prefix))
  }
}


################################################################################
# Manifest
################################################################################

generated <- sort(list.files(COR_DIR, full.names = TRUE, recursive = TRUE))
generated <- sub(paste0(COR_DIR, "/"), "", generated, fixed = TRUE)
writeLines(generated, file.path(COR_DIR, "AllUrothelium_correlation_files.txt"))

message(sprintf("\nOutputs saved to: %s", COR_DIR))
message(sprintf("Generated %d files across %d gene modes.",
                length(generated), length(GENE_MODES)))
