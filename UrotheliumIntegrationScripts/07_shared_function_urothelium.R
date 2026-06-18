################################################################################
# 07_shared_function_urothelium.R
#
# Goal: identify the pan-urothelial core gene program — genes and pathways
#       that are consistently expressed (shared function) across all urothelial
#       tissue contexts: kidney / bladder / ureter, in vivo + organoid.
#
# Input : output/AllUrothelium_gated_UMAP_nosnRNA.rds
#
# Outputs (output/AllUrothelium_shared_function/):
#   one_vs_rest/   — per-Category FindMarkers CSVs (upregulated vs. all others)
#   upset_upregulated_genes.pdf    — UpSet plot of gene-set overlaps
#   core_shared_genes.csv          — genes upregulated in >= MIN_SHARED_GROUPS
#   enrichment/    — GO BP + KEGG CSVs + dot-plots for core gene set
#   CoreUroScore_FeaturePlot.pdf   — UMAP coloured by core gene module score
#   CoreUroScore_ViolinPlot.pdf    — score distribution per tissue context
#   CoreUroScore_DotPlot.pdf       — average score per seurat_cluster
#
# Strategy
# ────────
#   1. One-vs-rest FindMarkers (Wilcoxon) per Category
#   2. Take genes significantly upregulated (log2FC > LFC_THRESH, adj.p < PADJ)
#      in each group vs. all others
#   3. Core genes = upregulated in >= MIN_SHARED_GROUPS groups
#   4. UpSet plot shows the full overlap landscape across all groups
#   5. GO BP + KEGG enrichment on core genes reveals shared biological functions
#   6. Module score (AddModuleScore) projects the core signature onto UMAP
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
})

options(future.globals.maxSize = 8 * 1024^3)

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR  <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/FinalUrotheliumCells"
SCR_DIR   <- file.path(BASE_DIR, "UrotheliumScripts")
OUT_DIR   <- file.path(SCR_DIR, "output")
SHR_DIR   <- file.path(OUT_DIR, "AllUrothelium_shared_function")
OVR_DIR   <- file.path(SHR_DIR, "one_vs_rest")
ENR_DIR   <- file.path(SHR_DIR, "enrichment")

for (d in c(SHR_DIR, OVR_DIR, ENR_DIR))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

args     <- commandArgs(trailingOnly = TRUE)
RDS_PATH <- if (length(args) >= 1) args[[1]] else
  file.path(OUT_DIR, "AllUrothelium_gated_UMAP.rds")

# ── Parameters ────────────────────────────────────────────────────────────────
LFC_THRESH        <- 0.5    # log2FC cutoff for one-vs-rest markers
PADJ_THRESH       <- 0.05   # adjusted p-value cutoff
MIN_PCT           <- 0.1    # min.pct for FindMarkers
MIN_SHARED_GROUPS <- 3      # gene must be significant in >= N groups → "core"
TOP_PATHWAY       <- 20     # top pathways shown in enrichment dot-plots
GO_PVAL_CUT       <- 0.05
KEGG_PVAL_CUT     <- 0.05

# ── Load ──────────────────────────────────────────────────────────────────────
message("Loading: ", RDS_PATH)
if (!file.exists(RDS_PATH)) stop("RDS not found: ", RDS_PATH)

so <- readRDS(RDS_PATH)
message(sprintf("  %d cells x %d genes", ncol(so), nrow(so)))

DefaultAssay(so) <- "RNA"
rna_layers <- Layers(so[["RNA"]])
if (!"data" %in% rna_layers) {
  message("  Joining layers to restore data layer ...")
  so <- JoinLayers(so)
}

message("  Categories: ", paste(sort(unique(so$Categories)), collapse = ", "))
message("  UMAP reduction available: ",
        paste(names(so@reductions), collapse = ", "))

# Pick the UMAP reduction name (handles umap_harmony or umap)
umap_key <- if ("umap_harmony" %in% names(so@reductions)) "umap_harmony" else "umap"


################################################################################
# Helper: GO BP + KEGG enrichment on a gene vector
################################################################################

run_enrichment <- function(genes, label, out_dir) {
  if (length(genes) < 5) {
    message(sprintf("  [enrich] %s: too few genes (%d) — skipping", label, length(genes)))
    return(invisible(NULL))
  }

  mapped <- tryCatch(
    bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db),
    error = function(e) NULL
  )
  if (is.null(mapped) || nrow(mapped) < 5) {
    message(sprintf("  [enrich] %s: Entrez mapping failed — skipping", label))
    return(invisible(NULL))
  }
  entrez <- mapped$ENTREZID

  # GO Biological Process
  go_res <- tryCatch(
    enrichGO(
      gene          = entrez,
      OrgDb         = org.Mm.eg.db,
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = GO_PVAL_CUT,
      readable      = TRUE
    ),
    error = function(e) { message("  GO failed: ", e$message); NULL }
  )

  if (!is.null(go_res) && nrow(go_res@result) > 0) {
    go_df <- as.data.frame(go_res)
    write.csv(go_df, file.path(out_dir, sprintf("%s_GO_BP.csv", label)),
              row.names = FALSE)
    p <- dotplot(go_res, showCategory = TOP_PATHWAY) +
      ggtitle(sprintf("%s — GO Biological Process", label))
    ggsave(file.path(out_dir, sprintf("%s_GO_BP_dotplot.pdf", label)),
           plot = p, width = 10, height = 9)
    message(sprintf("  [GO BP]  %s: %d terms", label, nrow(go_df)))
  }

  # KEGG
  kegg_res <- tryCatch(
    enrichKEGG(
      gene          = entrez,
      organism      = "mmu",
      pAdjustMethod = "BH",
      pvalueCutoff  = KEGG_PVAL_CUT
    ),
    error = function(e) { message("  KEGG failed: ", e$message); NULL }
  )

  if (!is.null(kegg_res) && nrow(kegg_res@result) > 0) {
    kegg_df <- as.data.frame(kegg_res)
    write.csv(kegg_df, file.path(out_dir, sprintf("%s_KEGG.csv", label)),
              row.names = FALSE)
    p <- dotplot(kegg_res, showCategory = TOP_PATHWAY) +
      ggtitle(sprintf("%s — KEGG Pathways", label))
    ggsave(file.path(out_dir, sprintf("%s_KEGG_dotplot.pdf", label)),
           plot = p, width = 10, height = 8)
    message(sprintf("  [KEGG]   %s: %d pathways", label, nrow(kegg_df)))
  }

  invisible(list(go = go_res, kegg = kegg_res))
}


################################################################################
# STEP 1: One-vs-rest FindMarkers per Category
################################################################################

message("\n========== STEP 1: One-vs-rest FindMarkers ==========")

categories <- sort(unique(as.character(so$Categories)))
message(sprintf("  Groups: %s", paste(categories, collapse = ", ")))
Idents(so) <- "Categories"

marker_sets <- setNames(
  lapply(categories, function(cat) {
    n <- sum(so$Categories == cat)
    if (n < 10) {
      message(sprintf("  Skipping %-20s: only %d cells", cat, n))
      return(NULL)
    }
    message(sprintf("  FindMarkers: %-20s (%d cells vs rest) ...", cat, n))
    tryCatch({
      m <- FindMarkers(
        so,
        ident.1         = cat,
        test.use        = "wilcox",
        logfc.threshold = LFC_THRESH,
        min.pct         = MIN_PCT,
        verbose         = FALSE
      )
      m$gene     <- rownames(m)
      m$category <- cat
      write.csv(m,
                file.path(OVR_DIR, sprintf("%s_vs_rest_markers.csv", cat)),
                row.names = FALSE)
      m
    }, error = function(e) {
      message(sprintf("  FindMarkers failed for %s: %s", cat, e$message))
      NULL
    })
  }),
  categories
)
marker_sets <- Filter(Negate(is.null), marker_sets)

# Significant upregulated genes per category
up_sets <- lapply(marker_sets, function(m) {
  m %>%
    filter(avg_log2FC > 0, p_val_adj < PADJ_THRESH) %>%
    pull(gene)
})

message("\n  Significant upregulated genes per group:")
for (cat in names(up_sets))
  message(sprintf("    %-22s: %d genes", cat, length(up_sets[[cat]])))


################################################################################
# STEP 2: UpSet plot — overlap landscape across all groups
################################################################################

message("\n========== STEP 2: UpSet plot ==========")

if (requireNamespace("UpSetR", quietly = TRUE) && length(up_sets) >= 2) {
  gene_universe <- sort(unique(unlist(up_sets)))
  upset_mat <- as.data.frame(
    lapply(up_sets, function(g) as.integer(gene_universe %in% g))
  )
  rownames(upset_mat) <- gene_universe

  pdf(file.path(SHR_DIR, "upset_upregulated_genes.pdf"), width = 14, height = 7)
  print(
    UpSetR::upset(
      upset_mat,
      sets             = rev(names(up_sets)),
      order.by         = "freq",
      decreasing       = TRUE,
      mainbar.y.label  = "Number of shared genes",
      sets.x.label     = "Upregulated genes per group",
      text.scale       = c(1.4, 1.2, 1.2, 1, 1.4, 1.2),
      point.size       = 3,
      line.size        = 1
    )
  )
  dev.off()
  message("  Saved: upset_upregulated_genes.pdf")
} else if (!requireNamespace("UpSetR", quietly = TRUE)) {
  message("  UpSetR not installed — skipping UpSet plot.")
  message("  Install with: install.packages('UpSetR')")
}


################################################################################
# STEP 3: Core genes — upregulated in >= MIN_SHARED_GROUPS groups
################################################################################

message(sprintf("\n========== STEP 3: Core genes (>= %d groups) ==========",
                MIN_SHARED_GROUPS))

gene_freq  <- sort(table(unlist(up_sets)), decreasing = TRUE)
core_genes <- names(gene_freq[gene_freq >= MIN_SHARED_GROUPS])

message(sprintf("  Core genes found: %d", length(core_genes)))

# Build summary table: which groups each core gene appears in
core_df <- data.frame(
  gene     = core_genes,
  n_groups = as.integer(gene_freq[core_genes]),
  groups   = sapply(core_genes, function(g)
    paste(names(up_sets)[sapply(up_sets, function(gs) g %in% gs)],
          collapse = ";")),
  stringsAsFactors = FALSE
) %>%
  arrange(desc(n_groups), gene)

write.csv(core_df,
          file.path(SHR_DIR, "core_shared_genes.csv"),
          row.names = FALSE)
message("  Saved: core_shared_genes.csv")

# Also save frequency table across all sharing levels
freq_summary <- as.data.frame(table(gene_freq), stringsAsFactors = FALSE)
colnames(freq_summary) <- c("n_groups_shared", "n_genes")
write.csv(freq_summary,
          file.path(SHR_DIR, "gene_sharing_frequency.csv"),
          row.names = FALSE)

# Bar plot: how many genes are shared at each threshold
p_bar <- ggplot(
  data.frame(n = as.integer(names(gene_freq)),
             gene = names(gene_freq)) %>%
    count(n, name = "count"),
  aes(x = factor(n), y = count)
) +
  geom_col(fill = "#2c7bb6", width = 0.6) +
  geom_text(aes(label = count), vjust = -0.4, size = 3.5) +
  labs(
    title = "Upregulated genes by number of groups sharing them",
    x     = "Number of tissue groups where gene is upregulated",
    y     = "Number of genes"
  ) +
  theme_bw(base_size = 12)
ggsave(file.path(SHR_DIR, "gene_sharing_barplot.pdf"),
       plot = p_bar, width = 7, height = 5)
message("  Saved: gene_sharing_barplot.pdf")


################################################################################
# STEP 4: GO BP + KEGG enrichment on core genes
################################################################################

message("\n========== STEP 4: Enrichment on core genes ==========")

run_enrichment(core_genes, "CoreUrothelium", ENR_DIR)

# Also run enrichment at the "all groups" level (genes in every group)
all_shared <- names(gene_freq[gene_freq == length(up_sets)])
if (length(all_shared) >= 5) {
  message(sprintf("  Genes shared in ALL %d groups: %d",
                  length(up_sets), length(all_shared)))
  write.csv(data.frame(gene = all_shared),
            file.path(SHR_DIR, "universal_shared_genes.csv"),
            row.names = FALSE)
  run_enrichment(all_shared, "UniversalUrothelium", ENR_DIR)
}


################################################################################
# STEP 5: Module score on UMAP + violin per tissue context
################################################################################

message("\n========== STEP 5: Core gene module score ==========")

genes_present <- intersect(core_genes, rownames(so))
message(sprintf("  Core genes in object: %d / %d", length(genes_present), length(core_genes)))

if (length(genes_present) >= 5) {
  so <- AddModuleScore(
    so,
    features = list(genes_present),
    name     = "CoreUroScore",
    seed     = 0
  )

  # FeaturePlot on UMAP
  p_feat <- FeaturePlot(
    so,
    features  = "CoreUroScore1",
    reduction = umap_key,
    raster    = TRUE,
    order     = TRUE,
    cols      = c("lightgrey", "#d73027")
  ) +
    ggtitle("Pan-urothelial core gene score") +
    theme(plot.title = element_text(size = 12))

  ggsave(file.path(SHR_DIR, "CoreUroScore_FeaturePlot.pdf"),
         plot = p_feat, width = 7, height = 6)
  message("  Saved: CoreUroScore_FeaturePlot.pdf")

  # Violin by Categories
  p_viol <- VlnPlot(
    so,
    features = "CoreUroScore1",
    group.by = "Categories",
    pt.size  = 0,
    cols     = scales::hue_pal()(length(unique(so$Categories)))
  ) +
    labs(
      title = "Pan-urothelial core gene score by tissue context",
      x     = NULL,
      y     = "Module score"
    ) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "none")

  ggsave(file.path(SHR_DIR, "CoreUroScore_ViolinPlot.pdf"),
         plot = p_viol, width = 9, height = 5)
  message("  Saved: CoreUroScore_ViolinPlot.pdf")

  # DimPlot by Categories with score overlay (split view)
  if ("seurat_clusters" %in% colnames(so@meta.data)) {
    p_dot <- DotPlot(
      so,
      features = "CoreUroScore1",
      group.by = "seurat_clusters"
    ) +
      RotatedAxis() +
      labs(title = "Core gene score — average per cluster") +
      theme(axis.text.x = element_text(size = 9))

    ggsave(file.path(SHR_DIR, "CoreUroScore_DotPlot_bycluster.pdf"),
           plot = p_dot, width = 10, height = 4)
    message("  Saved: CoreUroScore_DotPlot_bycluster.pdf")
  }

  # DimPlot split by Categories, coloured by score
  cats <- sort(unique(so$Categories))
  n_cat <- length(cats)
  p_split <- FeaturePlot(
    so,
    features  = "CoreUroScore1",
    reduction = umap_key,
    split.by  = "Categories",
    raster    = TRUE,
    order     = TRUE,
    cols      = c("lightgrey", "#d73027")
  )
  ggsave(file.path(SHR_DIR, "CoreUroScore_FeaturePlot_split.pdf"),
         plot = p_split,
         width = 5 * n_cat, height = 5)
  message("  Saved: CoreUroScore_FeaturePlot_split.pdf")

} else {
  message("  Too few core genes present in object — skipping module score.")
}


################################################################################
# STEP 6: Top core genes — dot-plot across all Categories
################################################################################

message("\n========== STEP 6: Core gene expression dot-plot ==========")

# Show the top genes shared in the most groups (max 30 for readability)
top_core <- head(core_df$gene[order(core_df$n_groups, decreasing = TRUE)], 30)
top_core <- intersect(top_core, rownames(so))

if (length(top_core) >= 2 && "Categories" %in% colnames(so@meta.data)) {
  p_core_dot <- DotPlot(
    so,
    features = top_core,
    group.by = "Categories",
    cols     = c("lightgrey", "#08519c"),
    dot.scale = 6
  ) +
    RotatedAxis() +
    labs(title = sprintf("Top %d core pan-urothelial genes", length(top_core))) +
    theme(axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 9))

  ggsave(file.path(SHR_DIR, "top_core_genes_DotPlot.pdf"),
         plot = p_core_dot,
         width = max(8, 2 + length(top_core) * 0.5),
         height = max(4, 2 + length(categories) * 0.4))
  message("  Saved: top_core_genes_DotPlot.pdf")
}


################################################################################
# Summary
################################################################################

message("\n========== Summary ==========")
message(sprintf("  Input cells              : %d", ncol(so)))
message(sprintf("  Groups analysed          : %d  (%s)",
                length(marker_sets), paste(names(marker_sets), collapse = ", ")))
message(sprintf("  Sharing threshold        : >= %d groups", MIN_SHARED_GROUPS))
message(sprintf("  Core genes               : %d", length(core_genes)))
if (length(all_shared) >= 5)
  message(sprintf("  Universal genes (all)    : %d", length(all_shared)))
message(sprintf("  Outputs                  : %s", SHR_DIR))

message("\nShared function analysis complete.")
