################################################################################
# 03_metadata_correlations_AllUrothelium.R
#
# Input : output/AllUrothelium_harmony_integrated.rds
#
# Outputs (output/AllUrothelium_metadata_correlations/):
#   1. Cramer's V association heatmap across categorical metadata columns
#   2. Cluster-composition Pearson correlation heatmaps per metadata column
#   3. Average-expression Pearson correlation heatmaps per metadata column
#   4. Supporting CSVs (matrices, summaries)
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(Matrix)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/FinalUrotheliumCells"
SCR_DIR  <- file.path(BASE_DIR, "UrotheliumScripts")
OUT_DIR  <- file.path(SCR_DIR, "output")
COR_DIR  <- file.path(OUT_DIR, "AllUrothelium_metadata_correlations")

args     <- commandArgs(trailingOnly = TRUE)
RDS_PATH <- if (length(args) >= 1) args[[1]] else
  file.path(OUT_DIR, "AllUrothelium_gated_UMAP.rds")
if (length(args) >= 2) COR_DIR <- args[[2]]

dir.create(COR_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Parameters 

────────────────────────────────────────────────────────────────
METADATA_COR_COLS   <- c("LabelClass","Finalgsm_id", "FinalConditionL1","FinalConditionL2", "sample_id", "technology", "paper", "tissue","Categories","Finalpaper")
CLUSTER_COL         <- "seurat_clusters"
MIN_CELLS_PER_LEVEL <- 20
MIN_EXPR_MEAN       <- 0.01   # filter genes with very low average expression
TOP_VAR_GENES       <- 2000   # use top variable genes for avg-expr correlation
DROP_CONDITIONS     <- c("E9To13.5Gestation", "E18_5_Kidney")

# ── Helpers ───────────────────────────────────────────────────────────────────

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

plot_correlation_heatmap <- function(mat, title, out_pdf, out_csv) {
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
    labs(title = title, x = NULL, y = NULL, fill = "Pearson r") +
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


################################################################################
# Load
################################################################################

message("Loading: ", RDS_PATH)
if (!file.exists(RDS_PATH))
  stop("Input RDS not found: ", RDS_PATH)

so   <- readRDS(RDS_PATH)
meta <- so@meta.data
message(sprintf("Loaded %s cells.", format(nrow(meta), big.mark = ",")))

# Drop developmental conditions
if ("condition" %in% colnames(meta)) {
  keep <- !(as.character(meta$condition) %in% DROP_CONDITIONS)
  if (any(!keep)) {
    so   <- so[, keep]
    meta <- so@meta.data
    message(sprintf("After dropping developmental conditions: %s cells.",
                    format(nrow(meta), big.mark = ",")))
  }
}

# Identify available metadata columns
metadata_cols <- intersect(METADATA_COR_COLS, colnames(meta))
if (length(metadata_cols) == 0)
  stop("None of the requested metadata columns found: ",
       paste(METADATA_COR_COLS, collapse = ", "))
if (!(CLUSTER_COL %in% colnames(meta)))
  stop("Cluster column not found: ", CLUSTER_COL)

# Metadata summary
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
# Average expression (computed before freeing the object)
################################################################################

message("Computing average expression per metadata group ...")

# Select genes for correlation: high-variance + minimum mean expression
rna_layers <- Layers(so[["RNA"]])
if (!"data" %in% rna_layers) so <- JoinLayers(so)

# Use variable features if available, otherwise all genes
hvg_genes <- VariableFeatures(so)
if (length(hvg_genes) == 0) {
  message("  No HVGs found — using all genes (may be slow).")
  hvg_genes <- rownames(so)
}
hvg_genes <- head(hvg_genes, TOP_VAR_GENES)

# Compute average log-normalised expression directly from the data layer.
# AverageExpression (Seurat v5) routes through PseudobulkExpression, which
# exponentiates the data layer expecting raw counts — yielding Inf for scVI
# log-norm values and mangling column names (underscores→dashes, numeric
# prefixes get 'g' prepended).  Direct rowMeans on the data matrix avoids all
# of that.
data_mat <- GetAssayData(so, assay = "RNA", layer = "data")
data_mat <- data_mat[intersect(hvg_genes, rownames(data_mat)), , drop = FALSE]

avg_expr_list <- setNames(
  lapply(metadata_cols, function(grp) {
    tryCatch({
      groups      <- as.character(so@meta.data[[grp]])
      valid_grps  <- names(which(table(groups[!is.na(groups)]) >= 1))
      if (length(valid_grps) < 2) return(NULL)
      avg_cols <- lapply(valid_grps, function(g) {
        idx <- which(!is.na(groups) & groups == g)
        Matrix::rowMeans(data_mat[, idx, drop = FALSE])
      })
      names(avg_cols) <- valid_grps
      do.call(cbind, avg_cols)
    }, error = function(e) {
      warning("avg-expr failed for ", grp, ": ", conditionMessage(e))
      NULL
    })
  }),
  metadata_cols
)
avg_expr_list <- Filter(Negate(is.null), avg_expr_list)

rm(so); gc()
message("  Average expression computed; Seurat object freed.")


################################################################################
# 1. Cramer's V association heatmap across metadata columns
################################################################################

if (length(metadata_cols) >= 2) {
  message("Computing Cramer's V association matrix ...")
  assoc_mat <- matrix(
    NA_real_,
    nrow     = length(metadata_cols),
    ncol     = length(metadata_cols),
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
# 2. Cluster-composition correlation heatmaps
################################################################################

message("Generating cluster-composition correlation heatmaps ...")

for (grp in metadata_cols) {
  meta_df <- meta[, c(CLUSTER_COL, grp), drop = FALSE]
  colnames(meta_df) <- c("cluster", "group")
  meta_df <- meta_df[!is.na(meta_df$cluster) & !is.na(meta_df$group), ,
                     drop = FALSE]
  meta_df$cluster <- as.character(meta_df$cluster)
  meta_df$group   <- as.character(meta_df$group)

  grp_counts  <- sort(table(meta_df$group), decreasing = TRUE)
  keep_groups <- names(grp_counts[grp_counts >= MIN_CELLS_PER_LEVEL])
  if (length(keep_groups) < 2) {
    warning("Skipping cluster-composition for ", grp,
            ": fewer than two levels with ≥ ", MIN_CELLS_PER_LEVEL, " cells.")
    next
  }

  meta_df <- meta_df[meta_df$group %in% keep_groups, , drop = FALSE]
  clust_levels <- sorted_levels(meta_df$cluster)
  grp_levels   <- keep_groups

  counts_mat <- table(
    cluster = factor(meta_df$cluster, levels = clust_levels),
    group   = factor(meta_df$group,   levels = grp_levels)
  )
  prop_mat <- prop.table(counts_mat, margin = 2)
  if (nrow(prop_mat) < 2) next

  cor_mat <- cor(as.matrix(prop_mat), method = "pearson",
                 use = "pairwise.complete.obs")

  prefix <- paste0("AllUrothelium_clustcomp_", safe_filename(grp))
  write.csv(as.data.frame.matrix(counts_mat),
            file.path(COR_DIR, paste0(prefix, "_counts.csv")), quote = FALSE)
  write.csv(as.data.frame.matrix(prop_mat),
            file.path(COR_DIR, paste0(prefix, "_proportions.csv")), quote = FALSE)
  plot_correlation_heatmap(
    cor_mat,
    sprintf("%s: cluster-composition correlation", grp),
    file.path(COR_DIR, paste0(prefix, "_correlation_heatmap.pdf")),
    file.path(COR_DIR, paste0(prefix, "_correlation_matrix.csv"))
  )
  message(sprintf("  Saved: %s_correlation_heatmap.pdf", prefix))
}


################################################################################
# 3. Average-expression correlation heatmaps
################################################################################

message("Generating average-expression correlation heatmaps ...")

for (grp in names(avg_expr_list)) {
  avg_mat <- avg_expr_list[[grp]]

  grp_counts  <- table(meta[[grp]][!is.na(meta[[grp]])])
  keep_groups <- names(grp_counts[grp_counts >= MIN_CELLS_PER_LEVEL])
  keep_groups <- intersect(keep_groups, colnames(avg_mat))
  if (length(keep_groups) < 2) {
    warning("Skipping avg-expr correlation for ", grp,
            ": fewer than two qualifying levels.")
    next
  }

  avg_mat <- avg_mat[, keep_groups, drop = FALSE]

  # Filter low-expressed genes
  gene_means <- rowMeans(as.matrix(avg_mat))
  avg_mat    <- avg_mat[gene_means >= MIN_EXPR_MEAN, , drop = FALSE]
  if (nrow(avg_mat) < 2) next

  cor_mat <- cor(as.matrix(avg_mat), method = "pearson",
                 use = "pairwise.complete.obs")

  prefix <- paste0("AllUrothelium_avgexpr_", safe_filename(grp))
  write.csv(avg_mat,
            file.path(COR_DIR, paste0(prefix, "_avg_expression.csv")),
            quote = FALSE)
  plot_correlation_heatmap(
    cor_mat,
    sprintf("%s: average-expression correlation", grp),
    file.path(COR_DIR, paste0(prefix, "_correlation_heatmap.pdf")),
    file.path(COR_DIR, paste0(prefix, "_correlation_matrix.csv"))
  )
  message(sprintf("  Saved: %s_correlation_heatmap.pdf", prefix))
}


################################################################################
# Manifest
################################################################################

generated <- sort(list.files(COR_DIR, full.names = FALSE))
writeLines(generated, file.path(COR_DIR, "AllUrothelium_correlation_files.txt"))

message(sprintf("\nMetadata correlation outputs saved to: %s", COR_DIR))
message(sprintf("Generated %d files.", length(generated)))











########### please ignore the following code, this is just for me to check the condition column in the metadata
# okay I need to filter this sample with development out: MultiOmicSpatialMouseKidney_SexSpecific_lifespan_Chen2025NatGenet

so@meta.data %>% filter(sample_id == "AKIGFPneg6mo") %>% head()
so@meta.data %>% filter(sample_id %in% c("AKIGFPneg4wk", "AKIGFPneg6mo", "AKIGFPpos4wk1", "AKIGFPpos4wk2", "AKIGFPpos4wk3", "AKIGFPpos6mo1", "AKIGFPpos6mo3", "AKIGFPpos6mo4","AMK13", "AMK14", "AMK16", "AMK17", "AMK18", "AMK27", "AMK29", "AMK33", "AMK34", "AMK35", "AMK36", "AMK37", "AMK38", "AMK39", "AMK4", "AMK40", "AMK5", "AMK6", "AMK8", "AMK9", "Cont4wk1", "Cont4wk2", "Cont6mo") ) %>% select(orig.ident, sample_id, condition, paper, tissue, technology, source_GEO,gsm_id, gse_id, Categories, scanvi_MKA_label, scanvi_Lake_label, condition_level1,  condition_level2, condition_level3) %>% filter(condition_level1 == "AKI") %>% pull(sample_id) %>% table() %>% sort(decreasing = TRUE)



# this is okay but I want to check the condition column too
so@meta.data %>% filter(sample_id %in% c("IRI12h1b1", "IRI12h1b2", "IRI12h2", "IRI12h3", "IRI14d1b1", "IRI14d1b2", "IRI14d2", "IRI14d3", "IRI2d1b1", "IRI2d1b2", "IRI2d2b1", "IRI2d2b2", "IRI2d3", "IRI4h1", "IRI4h2", "IRI4h3", "IRI6w1b1", "IRI6w1b2", "IRI6w2", "IRI6w3", "IRIsham1b1", "IRIsham1b2", "IRIsham2", "IRIsham3")) %>% select(orig.ident, sample_id, condition, paper, tissue, technology, source_GEO,gsm_id, gse_id, Categories, scanvi_MKA_label, scanvi_Lake_label, condition_level1,  condition_level2, condition_level3)  %>% pull(condition_level1) %>% table() %>% sort(decreasing = TRUE)


# this is okay but I want to check the condition column too
so@meta.data %>% filter(sample_id %in% c("IRI12h1b1", "IRI12h1b2", "IRI12h2", "IRI12h3", "IRI14d1b1", "IRI14d1b2", "IRI14d2", "IRI14d3", "IRI2d1b1", "IRI2d1b2", "IRI2d2b1", "IRI2d2b2", "IRI2d3", "IRI4h1", "IRI4h2", "IRI4h3", "IRI6w1b1", "IRI6w1b2", "IRI6w2", "IRI6w3", "IRIsham1b1", "IRIsham1b2", "IRIsham2", "IRIsham3")) %>% select(orig.ident, sample_id, condition, paper, tissue, technology, source_GEO,gsm_id, gse_id, Categories, scanvi_MKA_label, scanvi_Lake_label, condition_level1,  condition_level2, condition_level3)  %>% pull(condition) %>% table() %>% sort(decreasing = TRUE)

# need to change the condition from Shame to "Reference"
so@meta.data %>% filter(sample_id %in% c("KidneyHealthy3", "KidneyHealthy4", "KidneyHealthy5")) %>% select(orig.ident, sample_id, condition, paper, tissue, technology, source_GEO,gsm_id, gse_id, Categories, scanvi_MKA_label, scanvi_Lake_label, condition_level1,  condition_level2, condition_level3)  %>% pull(condition) %>% table() %>% sort(decreasing = TRUE)

# these samples are okay for the condition column
so@meta.data %>% filter(sample_id %in% c("KidneyUUO2", "KidneyUUO3", "KidneyUUO4", "KidneyUUO5", "KidneyUUO6")) %>% select(orig.ident, sample_id, condition, paper, tissue, technology, source_GEO,gsm_id, gse_id, Categories, scanvi_MKA_label, scanvi_Lake_label, condition_level1,  condition_level2, condition_level3)  %>% pull(condition) %>% table() %>% sort(decreasing = TRUE)

so@meta.data %>% filter(sample_id %in% c("D1", "D2")) %>% select(orig.ident, sample_id, condition, paper, tissue, technology, source_GEO,gsm_id, gse_id, Categories, scanvi_MKA_label, scanvi_Lake_label, condition_level1,  condition_level2, condition_level3)  %>% pull(condition) %>% table() %>% sort(decreasing = TRUE)

# this has some problems with the condition column
so@meta.data %>% filter(sample_id %in% c("UUO")) %>% select(orig.ident, sample_id, condition, paper, tissue, technology, source_GEO,gsm_id, gse_id, Categories, scanvi_MKA_label, scanvi_Lake_label, condition_level1,  condition_level2, condition_level3)  

# put into UUO sample to condition: GSE119531 , UUO_14days
so@meta.data %>% filter(sample_id %in% c("UUO"))%>% head()

# I want to put Urotehlium_Organoid seperate into several different conditions
so@meta.data %>% filter(condition %in% c("Urothelium_Organoid")) %>% select(orig.ident, sample_id, condition, paper, tissue, technology, source_GEO,gsm_id, gse_id, Categories, scanvi_MKA_label, scanvi_Lake_label, condition_level1,  condition_level2, condition_level3)  %>% head()

# these need to chagne the gsm_id 
so@meta.data %>% filter(condition %in% c("Healthy_Urothelium")) %>% select(orig.ident, sample_id, condition, paper, tissue, technology, source_GEO,gsm_id, gse_id, Categories, scanvi_MKA_label, scanvi_Lake_label, condition_level1,  condition_level2, condition_level3)  %>% head()

# paste(KidneyUrogenoid_Organoid and orig.ident)
so@meta.data %>% filter(condition %in% c("KidneyUrogenoid_Organoid")) %>% select(orig.ident, sample_id, condition, paper, tissue, technology, source_GEO,gsm_id, gse_id, Categories, scanvi_MKA_label, scanvi_Lake_label, condition_level1,  condition_level2, condition_level3)  %>% head()



# write a csv file for the sample_id and condition column to check the condition column
MetaDataNeedCorretion <- so@meta.data %>% select(orig.ident, sample_id, condition, paper, tissue, technology, source_GEO,gsm_id, gse_id, Categories, condition_level1,  condition_level2, condition_level3) %>% as.data.frame()

# remove the rownames
rownames(MetaDataNeedCorretion) <- NULL

MetaDataNeedCorretion2<- MetaDataNeedCorretion%>% unique()

write.csv(MetaDataNeedCorretion2, file.path(COR_DIR, "MetaDataNeedCorretion.csv"), quote = FALSE, row.names = FALSE)
