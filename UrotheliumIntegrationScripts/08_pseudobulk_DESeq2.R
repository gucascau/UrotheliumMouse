################################################################################
# 08_pseudobulk_DESeq2.R
#
# Pseudobulk DESeq2 differential expression between FinalConditionL2 groups,
# using sample_id as biological replicates, run separately within each
# FinalConditionL1.
#
# Strategy
# ────────
# For each FinalConditionL1 (e.g. InVivo, Organoid, …):
#   1. Subset cells belonging to that group.
#   2. Aggregate (sum) counts per sample_id → pseudobulk count matrix
#      (genes × samples).  Samples with fewer than MIN_CELLS_PER_SAMPLE cells
#      are excluded.
#   3. Run DESeq2 for all pairwise FinalConditionL2 comparisons within the
#      group (requires ≥ MIN_SAMPLES_PER_GROUP replicates on each side).
#   4. Save DEG CSVs + volcano PDFs.
#   5. GO BP + KEGG enrichment for each significant comparison.
#
# NOTE on raw counts: some samples in this atlas have log-normalised values
# in the counts slot (scVI-pipeline source).  The pseudobulk matrix is
# rounded to integers so DESeq2 can ingest it; DESeq2's size-factor
# normalisation then handles cross-sample scale differences.
#
# Input : output/AllUrothelium_gated_UMAP.rds   (or first CLI argument)
# Output: output/AllUrothelium_pseudobulk_DEG/
#           ├─ DEG/        — DESeq2 result CSVs
#           ├─ plots/      — volcano PDFs
#           └─ enrichment/ — GO/KEGG CSVs + dot-plots
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(Matrix)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/FinalUrotheliumCells"
SCR_DIR  <- file.path(BASE_DIR, "UrotheliumScripts")
OUT_DIR  <- file.path(SCR_DIR, "output")
PB_DIR   <- file.path(OUT_DIR, "AllUrothelium_pseudobulk_DEG")
DE_DIR   <- file.path(PB_DIR, "DEG")
PLT_DIR  <- file.path(PB_DIR, "plots")
ENR_DIR  <- file.path(PB_DIR, "enrichment")
for (d in c(DE_DIR, PLT_DIR, ENR_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

args   <- commandArgs(trailingOnly = TRUE)
RDS_IN <- if (length(args) >= 1) args[[1]] else
  file.path(OUT_DIR, "AllUrothelium_gated_UMAP.rds")

# ── Parameters ────────────────────────────────────────────────────────────────
L1_COL                <- "FinalConditionL1"  # outer grouping variable
L2_COL                <- "FinalConditionL2"  # comparison variable within each L1
SAMPLE_COL            <- "sample_id"         # biological replicate unit

MIN_CELLS_PER_SAMPLE  <- 10   # exclude samples with fewer cells than this
MIN_SAMPLES_PER_GROUP <- 2    # skip L2 groups with fewer samples than this
MIN_COUNT_FILTER      <- 10   # exclude genes with total pseudobulk counts < this

LFC_THRESH   <- 1.0    # |log2FC| threshold for volcano highlight
PADJ_THRESH  <- 0.05   # FDR threshold for DEG calling and enrichment
TOP_LABEL    <- 20     # number of gene labels on each volcano
GO_PVAL_CUT  <- 0.05
KEGG_PVAL_CUT <- 0.05
TOP_PATHWAY  <- 20


################################################################################
# Helper functions
################################################################################

safe_filename <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

# Aggregate counts per sample_id; returns a genes × samples dense integer matrix.
aggregate_pseudobulk <- function(cnt_mat, meta_sub, sample_col, min_cells) {
  samples <- unique(meta_sub[[sample_col]])
  pb_list <- lapply(samples, function(sid) {
    idx <- rownames(meta_sub)[meta_sub[[sample_col]] == sid]
    if (length(idx) < min_cells) {
      message(sprintf("    Skipping sample '%s': %d cells (< %d required)",
                      sid, length(idx), min_cells))
      return(NULL)
    }
    Matrix::rowSums(cnt_mat[, idx, drop = FALSE])
  })
  names(pb_list) <- samples
  pb_list <- Filter(Negate(is.null), pb_list)
  if (length(pb_list) == 0) return(NULL)
  mat <- do.call(cbind, pb_list)   # genes × samples (dense)
  mat <- round(mat)
  storage.mode(mat) <- "integer"
  mat
}

# Volcano plot saved as PDF.
plot_volcano <- function(res_df, title, out_pdf) {
  df <- res_df %>%
    mutate(
      gene     = rownames(.),
      sig      = !is.na(padj) & padj <= PADJ_THRESH & abs(log2FoldChange) >= LFC_THRESH,
      neg_logp = -log10(pmax(padj, 1e-300))
    )
  top_genes <- df %>% filter(sig) %>% arrange(padj) %>%
    slice_head(n = TOP_LABEL) %>% pull(gene)

  p <- ggplot(df, aes(x = log2FoldChange, y = neg_logp)) +
    geom_point(aes(color = sig), size = 0.8, alpha = 0.6, shape = 16) +
    scale_color_manual(values = c("FALSE" = "grey65", "TRUE" = "firebrick3"),
                       guide  = "none") +
    geom_vline(xintercept = c(-LFC_THRESH, LFC_THRESH),
               linetype = "dashed", color = "navy", linewidth = 0.5) +
    geom_hline(yintercept = -log10(PADJ_THRESH),
               linetype = "dashed", color = "navy", linewidth = 0.5) +
    ggrepel::geom_text_repel(
      data        = df %>% filter(gene %in% top_genes),
      aes(label   = gene),
      size        = 2.5,
      max.overlaps = 25,
      segment.size = 0.25
    ) +
    labs(
      title = title,
      x     = expression(log[2]~Fold~Change),
      y     = expression(-log[10](adjusted~italic(p)))
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(size = 10))

  ggsave(out_pdf, plot = p, width = 7, height = 6)
  invisible(p)
}

# GO BP + KEGG enrichment on significant DEGs.
run_enrichment <- function(res_df, tag, out_dir) {
  sig_genes <- rownames(res_df)[
    !is.na(res_df$padj) & res_df$padj <= PADJ_THRESH &
      abs(res_df$log2FoldChange) >= LFC_THRESH
  ]
  if (length(sig_genes) < 5) {
    message(sprintf("    Too few DEGs for enrichment (%d) — skipping", length(sig_genes)))
    return(invisible(NULL))
  }

  eg <- tryCatch(
    bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db),
    error = function(e) NULL
  )
  if (is.null(eg) || nrow(eg) == 0) {
    message("    bitr() returned no Entrez IDs — skipping enrichment")
    return(invisible(NULL))
  }
  entrez <- eg$ENTREZID

  # GO Biological Process
  go <- tryCatch(
    enrichGO(gene = entrez, OrgDb = org.Mm.eg.db, ont = "BP",
             readable = TRUE, pvalueCutoff = GO_PVAL_CUT, qvalueCutoff = 0.2),
    error = function(e) NULL
  )
  if (!is.null(go) && nrow(as.data.frame(go)) > 0) {
    write.csv(as.data.frame(go),
              file.path(out_dir, paste0("GO_BP_", tag, ".csv")), quote = FALSE)
    p <- dotplot(go, showCategory = TOP_PATHWAY) +
      ggtitle(paste("GO BP:", gsub("_", " ", tag)))
    ggsave(file.path(out_dir, paste0("GO_BP_", tag, ".pdf")), p, width = 8, height = 7)
    message(sprintf("    GO BP: %d terms", nrow(as.data.frame(go))))
  }

  # KEGG
  kegg <- tryCatch(
    enrichKEGG(gene = entrez, organism = "mmu",
               pvalueCutoff = KEGG_PVAL_CUT, qvalueCutoff = 0.2),
    error = function(e) NULL
  )
  if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
    write.csv(as.data.frame(kegg),
              file.path(out_dir, paste0("KEGG_", tag, ".csv")), quote = FALSE)
    p <- dotplot(kegg, showCategory = TOP_PATHWAY) +
      ggtitle(paste("KEGG:", gsub("_", " ", tag)))
    ggsave(file.path(out_dir, paste0("KEGG_", tag, ".pdf")), p, width = 8, height = 7)
    message(sprintf("    KEGG: %d pathways", nrow(as.data.frame(kegg))))
  }
}


################################################################################
# STEP 1: Load Seurat object
################################################################################

message("Loading: ", RDS_IN)
if (!file.exists(RDS_IN)) stop("RDS not found: ", RDS_IN)
so <- readRDS(RDS_IN)
message(sprintf("  %d cells × %d genes", ncol(so), nrow(so)))
DefaultAssay(so) <- "RNA"

# Join split layers (Seurat v5)
if (length(Layers(so[["RNA"]])) > 2) {
  message("  Joining RNA layers ...")
  so <- JoinLayers(so)
}

message(sprintf("  %s levels: %s",
                L1_COL, paste(sort(unique(so@meta.data[[L1_COL]])), collapse = ", ")))
message(sprintf("  %s levels: %s",
                L2_COL, paste(sort(unique(so@meta.data[[L2_COL]])), collapse = ", ")))


################################################################################
# STEP 2: Extract count matrix and cell metadata
################################################################################

message("Extracting count matrix ...")
cnt  <- GetAssayData(so, layer = "counts")   # genes × cells (sparse)
meta <- so@meta.data


################################################################################
# STEP 3: Pseudobulk DESeq2 per FinalConditionL1
################################################################################

l1_levels <- sort(unique(na.omit(meta[[L1_COL]])))
message(sprintf("\nFinalConditionL1 groups to process: %s",
                paste(l1_levels, collapse = ", ")))

for (l1 in l1_levels) {

  message(sprintf("\n%s\nProcessing FinalConditionL1 = '%s'", strrep("=", 60), l1))

  # ── Subset to this L1 ──────────────────────────────────────────────────────
  idx_l1  <- rownames(meta)[!is.na(meta[[L1_COL]]) & meta[[L1_COL]] == l1]
  meta_l1 <- meta[idx_l1, ]
  cnt_l1  <- cnt[, idx_l1, drop = FALSE]

  l2_levels <- sort(unique(na.omit(meta_l1[[L2_COL]])))
  message(sprintf("  FinalConditionL2 levels: %s", paste(l2_levels, collapse = ", ")))

  if (length(l2_levels) < 2) {
    message("  Skipping: fewer than 2 FinalConditionL2 groups.")
    next
  }

  # ── Pseudobulk aggregation ─────────────────────────────────────────────────
  pb_mat <- aggregate_pseudobulk(cnt_l1, meta_l1, SAMPLE_COL, MIN_CELLS_PER_SAMPLE)
  if (is.null(pb_mat)) {
    message("  Skipping: no samples passed cell-count QC.")
    next
  }

  # Filter low-count genes
  keep_genes <- rowSums(pb_mat) >= MIN_COUNT_FILTER
  pb_mat     <- pb_mat[keep_genes, , drop = FALSE]
  message(sprintf("  Genes after count filter: %d", nrow(pb_mat)))

  # Sample-level metadata (one row per sample_id)
  sample_meta <- meta_l1[
    match(colnames(pb_mat), meta_l1[[SAMPLE_COL]]), c(SAMPLE_COL, L2_COL)
  ]
  rownames(sample_meta) <- sample_meta[[SAMPLE_COL]]

  # ── All pairwise L2 comparisons ────────────────────────────────────────────
  pairs <- combn(l2_levels, 2, simplify = FALSE)

  for (pair in pairs) {
    grp_a <- pair[1]; grp_b <- pair[2]
    message(sprintf("\n  Comparing: '%s' vs '%s'  (positive LFC = higher in %s)",
                    grp_b, grp_a, grp_b))

    # Subset samples to this pair
    keep_samp  <- sample_meta[[L2_COL]] %in% c(grp_a, grp_b)
    pair_meta  <- sample_meta[keep_samp, , drop = FALSE]
    pair_mat   <- pb_mat[, rownames(pair_meta), drop = FALSE]

    # Check minimum replicates
    n_per_grp <- table(pair_meta[[L2_COL]])
    if (any(n_per_grp < MIN_SAMPLES_PER_GROUP)) {
      message(sprintf("    Skipping: insufficient replicates — %s",
                      paste(paste0(names(n_per_grp), "=", n_per_grp), collapse = ", ")))
      next
    }
    message(sprintf("    Samples per group: %s",
                    paste(paste0(names(n_per_grp), "=", n_per_grp), collapse = ", ")))

    pair_meta[[L2_COL]] <- factor(pair_meta[[L2_COL]], levels = c(grp_a, grp_b))

    # ── DESeq2 ──────────────────────────────────────────────────────────────
    res_df <- tryCatch({
      dds <- DESeqDataSetFromMatrix(
        countData = pair_mat,
        colData   = pair_meta,
        design    = as.formula(paste("~", L2_COL))
      )
      dds <- DESeq(dds, quiet = TRUE)
      res <- results(dds,
                     contrast = c(L2_COL, grp_b, grp_a),
                     alpha    = PADJ_THRESH)
      as.data.frame(res) %>%
        filter(!is.na(padj)) %>%
        arrange(padj)
    }, error = function(e) {
      message("    DESeq2 failed: ", conditionMessage(e))
      NULL
    })

    if (is.null(res_df)) next

    n_sig <- sum(res_df$padj <= PADJ_THRESH, na.rm = TRUE)
    message(sprintf("    Significant DEGs (padj ≤ %.2f, |LFC| ≥ %.1f): %d",
                    PADJ_THRESH, LFC_THRESH, n_sig))

    tag      <- sprintf("%s_%s_vs_%s",
                        safe_filename(l1), safe_filename(grp_b), safe_filename(grp_a))
    csv_path <- file.path(DE_DIR,  paste0("DEG_", tag, ".csv"))
    pdf_path <- file.path(PLT_DIR, paste0("volcano_", tag, ".pdf"))

    write.csv(res_df, csv_path, quote = FALSE)
    message("    Saved: ", csv_path)

    plot_volcano(res_df,
                 title   = sprintf("%s: %s vs %s", l1, grp_b, grp_a),
                 out_pdf = pdf_path)
    message("    Saved: ", pdf_path)

    run_enrichment(res_df, tag, ENR_DIR)
  }
}

message("\nDone.")
