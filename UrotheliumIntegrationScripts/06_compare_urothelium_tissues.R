################################################################################
# 06_compare_urothelium_tissues.R
#
# Differential expression + pathway enrichment comparing:
#
#   COMPARISON 1 — In vivo tissue origin
#     kidneyMouse  vs  bladderMouse
#     → What genes and pathways are specific to kidney urothelium in vivo?
#
#   COMPARISON 2 — Organoid tissue origin (3-way pairwise)
#     kidneyOrganoid  vs  bladderOrganoid
#     kidneyOrganoid  vs  ureterOrganoid
#     bladderOrganoid vs  ureterOrganoid
#
# For each comparison:
#   1. FindMarkers (Wilcoxon, log2FC ≥ 0.5, adj.p ≤ 0.05)
#   2. Volcano plot
#   3. Top-50 DE gene heatmap
#   4. GO Biological Process enrichment  (clusterProfiler)
#   5. KEGG pathway enrichment           (clusterProfiler)
#   6. Pathway dot-plot
#   7. PROGENy signaling pathway activity (14 pathways, Wilcoxon + violin + heatmap)
#
# Input : output/AllUrothelium_gated_UMAP.rds
# Output: output/AllUrothelium_tissue_comparison/
#           ├─ DE/          — marker CSVs
#           ├─ plots/       — volcano, DE heatmap, PROGENy violin + heatmap PDFs
#           └─ enrichment/  — GO / KEGG CSVs + dot-plots, PROGENy stats CSVs
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(clusterProfiler)
  library(org.Mm.eg.db)     # mouse gene annotation
  library(enrichplot)
  library(progeny)           # signaling pathway activity scoring
})

options(future.globals.maxSize = 8 * 1024^3)

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/FinalUrotheliumCells"
SCR_DIR  <- file.path(BASE_DIR, "UrotheliumScripts")
OUT_DIR  <- file.path(SCR_DIR, "output")
CMP_DIR  <- file.path(OUT_DIR, "AllUrothelium_tissue_comparison")
DE_DIR   <- file.path(CMP_DIR, "DE")
PLT_DIR  <- file.path(CMP_DIR, "plots")
ENR_DIR  <- file.path(CMP_DIR, "enrichment")
for (d in c(DE_DIR, PLT_DIR, ENR_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

args     <- commandArgs(trailingOnly = TRUE)
RDS_IN   <- if (length(args) >= 1) args[[1]] else
  file.path(OUT_DIR, "AllUrothelium_gated_UMAP.rds")

# ── Parameters ────────────────────────────────────────────────────────────────
LFC_THRESH   <- 0.5       # log2FC cutoff for FindMarkers
PADJ_THRESH  <- 0.05      # adjusted p-value cutoff
MIN_PCT      <- 0.1       # min.pct for FindMarkers
TOP_HEATMAP  <- 50        # top DE genes per direction for heatmap
TOP_PATHWAY  <- 20        # top pathways to display
GO_PVAL_CUT  <- 0.05
KEGG_PVAL_CUT <- 0.05

# ── Load ──────────────────────────────────────────────────────────────────────
message("Loading: ", RDS_IN)
if (!file.exists(RDS_IN)) stop("RDS not found: ", RDS_IN)
so <- readRDS(RDS_IN)
message(sprintf("  %d cells × %d genes", ncol(so), nrow(so)))
message("  Categories: ", paste(sort(unique(so$Categories)), collapse = ", "))
DefaultAssay(so) <- "RNA"


################################################################################
# Helper functions
################################################################################

# Convert gene symbols to Entrez IDs (needed for KEGG)
symbols_to_entrez <- function(genes) {
  mapped <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Mm.eg.db)
  mapped$ENTREZID
}

# Run GO BP + KEGG enrichment on a vector of gene symbols
run_enrichment <- function(genes_up, genes_dn, label, out_dir) {
  results <- list()

  for (direction in c("up", "dn")) {
    genes <- if (direction == "up") genes_up else genes_dn
    tag   <- if (direction == "up") "UP_in_ident1" else "DOWN_in_ident1"

    if (length(genes) < 5) {
      message(sprintf("  [%s %s] too few genes (%d) — skipping enrichment",
                      label, tag, length(genes)))
      next
    }

    entrez <- tryCatch(symbols_to_entrez(genes), error = function(e) NULL)
    if (is.null(entrez) || length(entrez) < 5) {
      message(sprintf("  [%s %s] Entrez mapping failed — skipping", label, tag))
      next
    }

    # GO Biological Process
    go_res <- tryCatch(
      enrichGO(gene          = entrez,
               OrgDb         = org.Mm.eg.db,
               ont           = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff  = GO_PVAL_CUT,
               readable      = TRUE),
      error = function(e) { message("  GO failed: ", e$message); NULL }
    )

    if (!is.null(go_res) && nrow(go_res@result) > 0) {
      go_df <- as.data.frame(go_res)
      write.csv(go_df, file.path(out_dir, sprintf("%s_%s_GO_BP.csv", label, tag)),
                row.names = FALSE)

      p <- dotplot(go_res, showCategory = TOP_PATHWAY) +
        ggtitle(sprintf("%s — %s — GO BP", label, tag))
      ggsave(file.path(out_dir, sprintf("%s_%s_GO_BP_dotplot.pdf", label, tag)),
             plot = p, width = 10, height = 8)
      results[[paste(direction, "go")]] <- go_df
      message(sprintf("  GO BP [%s %s]: %d terms", label, tag, nrow(go_df)))
    }

    # KEGG
    kegg_res <- tryCatch(
      enrichKEGG(gene          = entrez,
                 organism      = "mmu",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = KEGG_PVAL_CUT),
      error = function(e) { message("  KEGG failed (network?): ", e$message); NULL }
    )

    if (!is.null(kegg_res) && nrow(kegg_res@result) > 0) {
      kegg_df <- as.data.frame(kegg_res)
      write.csv(kegg_df, file.path(out_dir, sprintf("%s_%s_KEGG.csv", label, tag)),
                row.names = FALSE)

      p <- dotplot(kegg_res, showCategory = TOP_PATHWAY) +
        ggtitle(sprintf("%s — %s — KEGG", label, tag))
      ggsave(file.path(out_dir, sprintf("%s_%s_KEGG_dotplot.pdf", label, tag)),
             plot = p, width = 10, height = 7)
      results[[paste(direction, "kegg")]] <- kegg_df
      message(sprintf("  KEGG [%s %s]: %d pathways", label, tag, nrow(kegg_df)))
    }
  }
  invisible(results)
}


# Volcano plot
plot_volcano <- function(markers_df, label, lfc = LFC_THRESH, padj = PADJ_THRESH) {
  df <- markers_df %>%
    mutate(
      sig = case_when(
        avg_log2FC >= lfc  & p_val_adj < padj ~ "UP",
        avg_log2FC <= -lfc & p_val_adj < padj ~ "DOWN",
        TRUE ~ "NS"
      ),
      neg_log10_p = -log10(pmax(p_val_adj, 1e-300))
    )

  top_up   <- df %>% filter(sig == "UP")   %>% slice_max(avg_log2FC,       n = 15)
  top_dn   <- df %>% filter(sig == "DOWN") %>% slice_min(avg_log2FC,       n = 15)
  top_sig  <- df %>% filter(sig != "NS")   %>% slice_min(p_val_adj,        n = 20)
  label_df <- bind_rows(top_up, top_dn, top_sig) %>% distinct()

  n_up <- sum(df$sig == "UP");  n_dn <- sum(df$sig == "DOWN")

  ggplot(df, aes(avg_log2FC, neg_log10_p, colour = sig)) +
    geom_point(alpha = 0.5, size = 0.8) +
    scale_colour_manual(values = c(UP = "#d73027", DOWN = "#4575b4", NS = "grey70")) +
    geom_vline(xintercept = c(-lfc, lfc), linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = -log10(padj),  linetype = "dashed", colour = "grey40") +
    ggrepel::geom_text_repel(
      data = label_df, aes(label = gene),
      size = 2.5, max.overlaps = 20, colour = "black"
    ) +
    labs(
      title    = label,
      subtitle = sprintf("UP: %d  |  DOWN: %d", n_up, n_dn),
      x        = "avg log2FC  (ident1 vs ident2)",
      y        = "-log10(adj p-value)",
      colour   = NULL
    ) +
    theme_bw(base_size = 11)
}


# Top-N heatmap of DE genes
plot_heatmap <- function(so_sub, markers_df, label, n = TOP_HEATMAP) {
  top_genes <- markers_df %>%
    filter(p_val_adj < PADJ_THRESH) %>%
    {bind_rows(
      slice_max(., avg_log2FC, n = n),
      slice_min(., avg_log2FC, n = n)
    )} %>%
    pull(gene) %>% unique() %>%
    intersect(rownames(so_sub))

  if (length(top_genes) < 2) {
    message("  Too few genes for heatmap: ", label)
    return(invisible(NULL))
  }

  mat <- AverageExpression(so_sub, features = top_genes,
                           group.by = "plot_group",
                           assays = "RNA", layer = "data")$RNA
  mat <- t(scale(t(as.matrix(mat))))
  mat[is.nan(mat)] <- 0

  ann_col <- data.frame(Group = colnames(mat), row.names = colnames(mat))

  pheatmap(
    mat,
    annotation_col  = ann_col,
    cluster_rows    = TRUE,
    cluster_cols    = FALSE,
    show_rownames   = TRUE,
    show_colnames   = TRUE,
    fontsize_row    = 7,
    main            = label,
    filename        = file.path(PLT_DIR, sprintf("%s_heatmap.pdf", label)),
    width = 8, height = max(6, length(top_genes) * 0.18 + 2)
  )
  message(sprintf("  Heatmap saved: %s", label))
}


# PROGENy pathway activity comparison
run_progeny_comparison <- function(so_sub, ident1, ident2, label) {
  message(sprintf("  [PROGENy] %s", label))

  expr_mat <- GetAssayData(so_sub, assay = "RNA", layer = "data")

  scores <- tryCatch(
    progeny(
      expr     = as.matrix(expr_mat),
      scale    = TRUE,
      organism = "Mouse",
      top      = 500,
      perm     = 1,
      verbose  = FALSE
    ),
    error = function(e) {
      message("  [PROGENy] failed: ", e$message)
      return(NULL)
    }
  )
  if (is.null(scores)) return(invisible(NULL))

  # scores is cells × pathways — add group label
  scores_df <- as.data.frame(scores)
  scores_df$plot_group <- so_sub$plot_group

  pathways <- setdiff(colnames(scores_df), "plot_group")

  # Wilcoxon test per pathway
  stats <- lapply(pathways, function(pw) {
    g1 <- scores_df[scores_df$plot_group == ident1, pw]
    g2 <- scores_df[scores_df$plot_group == ident2, pw]
    wt <- wilcox.test(g1, g2, exact = FALSE)
    data.frame(
      pathway      = pw,
      median_ident1 = median(g1),
      median_ident2 = median(g2),
      delta_median  = median(g1) - median(g2),
      p_value       = wt$p.value,
      stringsAsFactors = FALSE
    )
  })
  stats_df <- bind_rows(stats) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    arrange(p_adj)

  write.csv(stats_df,
            file.path(ENR_DIR, sprintf("%s_PROGENy_stats.csv", label)),
            row.names = FALSE)
  message(sprintf("  [PROGENy] sig pathways (adj.p < 0.05): %d",
                  sum(stats_df$p_adj < 0.05, na.rm = TRUE)))

  # Violin plot — one panel per pathway
  plot_df <- scores_df %>%
    pivot_longer(all_of(pathways), names_to = "pathway", values_to = "score") %>%
    mutate(plot_group = factor(plot_group, levels = c(ident1, ident2)))

  sig_paths <- stats_df$pathway[stats_df$p_adj < 0.05]

  p_violin <- ggplot(plot_df, aes(plot_group, score, fill = plot_group)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.8) +
    geom_boxplot(width = 0.12, outlier.size = 0.3, fill = "white") +
    facet_wrap(~ pathway, scales = "free_y", ncol = 5) +
    scale_fill_manual(values = c("#d73027", "#4575b4")) +
    labs(title = sprintf("PROGENy pathway activity — %s", label),
         x = NULL, y = "Activity score", fill = NULL) +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          strip.text  = element_text(size = 8))

  ggsave(file.path(PLT_DIR, sprintf("%s_PROGENy_violin.pdf", label)),
         plot = p_violin, width = 16, height = 10)

  # Heatmap of average pathway activity per group
  avg_scores <- scores_df %>%
    group_by(plot_group) %>%
    summarise(across(all_of(pathways), mean), .groups = "drop") %>%
    tibble::column_to_rownames("plot_group") %>%
    t()                              # pathways × groups

  ann_col <- data.frame(Group = colnames(avg_scores),
                        row.names = colnames(avg_scores))
  row_ann <- data.frame(
    Significant = ifelse(rownames(avg_scores) %in% sig_paths, "yes", "no"),
    row.names   = rownames(avg_scores)
  )

  pheatmap(
    avg_scores,
    annotation_col  = ann_col,
    annotation_row  = row_ann,
    annotation_colors = list(Significant = c(yes = "#d73027", no = "grey80")),
    cluster_cols    = FALSE,
    cluster_rows    = TRUE,
    scale           = "row",
    fontsize_row    = 9,
    main            = sprintf("PROGENy — %s", label),
    filename        = file.path(PLT_DIR,
                                sprintf("%s_PROGENy_heatmap.pdf", label)),
    width = 6, height = 7
  )

  message(sprintf("  [PROGENy] plots saved for %s", label))
  invisible(stats_df)
}


# Run a full comparison: DE + volcano + heatmap + enrichment + PROGENy
run_comparison <- function(so_sub, ident1, ident2, label) {
  message(sprintf("\n========== %s ==========", label))

  n1 <- sum(so_sub$plot_group == ident1)
  n2 <- sum(so_sub$plot_group == ident2)
  message(sprintf("  %s: %d cells  |  %s: %d cells", ident1, n1, ident2, n2))

  if (n1 < 10 || n2 < 10) {
    message("  Skipping — too few cells (< 10) in one group")
    return(invisible(NULL))
  }

  Idents(so_sub) <- "plot_group"

  markers <- FindMarkers(
    so_sub,
    ident.1        = ident1,
    ident.2        = ident2,
    test.use       = "wilcox",
    logfc.threshold = LFC_THRESH,
    min.pct        = MIN_PCT,
    verbose        = FALSE
  )
  markers$gene <- rownames(markers)
  message(sprintf("  DE genes (|log2FC| ≥ %.1f, adj.p < %.2f): %d up / %d down",
                  LFC_THRESH, PADJ_THRESH,
                  sum(markers$avg_log2FC > 0  & markers$p_val_adj < PADJ_THRESH),
                  sum(markers$avg_log2FC < 0  & markers$p_val_adj < PADJ_THRESH)))

  # Save DE table
  write.csv(markers, file.path(DE_DIR, sprintf("%s_markers.csv", label)),
            row.names = FALSE)

  # Volcano
  p_vol <- plot_volcano(markers, label)
  ggsave(file.path(PLT_DIR, sprintf("%s_volcano.pdf", label)),
         plot = p_vol, width = 8, height = 6)

  # Heatmap
  plot_heatmap(so_sub, markers, label)

  # Enrichment
  genes_up <- markers %>%
    filter(avg_log2FC > 0 & p_val_adj < PADJ_THRESH) %>%
    pull(gene)
  genes_dn <- markers %>%
    filter(avg_log2FC < 0 & p_val_adj < PADJ_THRESH) %>%
    pull(gene)

  run_enrichment(genes_up, genes_dn, label, ENR_DIR)

  # PROGENy pathway activity
  run_progeny_comparison(so_sub, ident1, ident2, label)

  invisible(markers)
}


################################################################################
# COMPARISON 1: In vivo — kidneyMouse vs bladderMouse
################################################################################

message("\n\n########## COMPARISON 1: In vivo kidney vs bladder ##########")

invivo_cats <- c("kidneyMouse", "bladderMouse")

if (all(invivo_cats %in% so$Categories)) {
  so_invivo <- so[, so$Categories %in% invivo_cats]
  so_invivo$plot_group <- so_invivo$Categories
  message(sprintf("  In vivo subset: %d cells", ncol(so_invivo)))

  markers_invivo <- run_comparison(
    so_invivo,
    ident1 = "kidneyMouse",
    ident2 = "bladderMouse",
    label  = "KidneyInvivo_vs_BladderInvivo"
  )

  # Additional: stratify by condition to check if markers hold across disease states
  if (!is.null(markers_invivo)) {
    conditions <- so_invivo$FinalConditionL1 %>% unique() %>% na.omit()
    message(sprintf("  Conditions present: %s", paste(conditions, collapse = ", ")))

    # Healthy/reference only
    ref_conds <- grep("Healthy|Ref|Sham|Normal|WT|Control|Ctrl",
                      conditions, value = TRUE, ignore.case = TRUE)
    if (length(ref_conds) > 0) {
      so_ref <- so_invivo[, so_invivo$FinalConditionL1 %in% ref_conds]
      message(sprintf("  Healthy-only subset: %d cells", ncol(so_ref)))
      if (ncol(so_ref) > 50) {
        run_comparison(
          so_ref,
          ident1 = "kidneyMouse",
          ident2 = "bladderMouse",
          label  = "KidneyInvivo_vs_BladderInvivo_HealthyOnly"
        )
      }
    }
  }
} else {
  missing <- setdiff(invivo_cats, unique(so$Categories))
  message("  WARNING: missing Categories: ", paste(missing, collapse = ", "))
}


################################################################################
# COMPARISON 2: Organoid — pairwise (kidney vs bladder vs ureter)
################################################################################

message("\n\n########## COMPARISON 2: Organoid pairwise ##########")

organoid_cats <- c("kidneyOrganoid", "bladderOrganoid", "ureterOrganoid")
present_org   <- intersect(organoid_cats, unique(so$Categories))
message("  Organoid groups found: ", paste(present_org, collapse = ", "))

if (length(present_org) >= 2) {
  so_org <- so[, so$Categories %in% present_org]
  so_org$plot_group <- so_org$Categories
  message(sprintf("  Organoid subset: %d cells", ncol(so_org)))
  print(table(so_org$Categories))

  # All pairwise comparisons
  pairs <- combn(present_org, 2, simplify = FALSE)
  for (pair in pairs) {
    run_comparison(
      so_org,
      ident1 = pair[1],
      ident2 = pair[2],
      label  = paste0(pair[1], "_vs_", pair[2])
    )
  }
} else {
  message("  WARNING: fewer than 2 organoid groups found — skipping organoid comparison")
}


################################################################################
# COMPARISON 3: In vivo vs Organoid within kidney (optional cross-context check)
################################################################################

message("\n\n########## COMPARISON 3: Kidney invivo vs Kidney organoid ##########")

kidney_cats <- c("kidneyMouse", "kidneyOrganoid")
if (all(kidney_cats %in% so$Categories)) {
  so_kidney <- so[, so$Categories %in% kidney_cats]
  so_kidney$plot_group <- so_kidney$Categories

  run_comparison(
    so_kidney,
    ident1 = "kidneyMouse",
    ident2 = "kidneyOrganoid",
    label  = "KidneyInvivo_vs_KidneyOrganoid"
  )
} else {
  message("  Skipping — one or both kidney groups missing")
}


################################################################################
# COMPARISON 4: Healthy-only — kidney vs bladder in vivo
# COMPARISON 5: Healthy-only — organoid pairwise (kidney / bladder / ureter)
#
# "Healthy" is defined as cells whose FinalConditionL1 or condition column
# matches any of: Healthy, Normal, WT, Ref, Sham, Control, Ctrl
################################################################################

healthy_pattern <- "Healthy|Normal|WT|Ref|Sham|Control|Ctrl"

# Identify healthy cells using whichever condition column is available
cond_col <- if ("FinalConditionL1" %in% colnames(so@meta.data)) "FinalConditionL1" else "condition"
message(sprintf("\nUsing '%s' to define healthy cells (pattern: %s)", cond_col, healthy_pattern))

is_healthy <- grepl(healthy_pattern, so@meta.data[[cond_col]], ignore.case = TRUE)
message(sprintf("  Healthy cells: %d / %d total", sum(is_healthy), ncol(so)))
print(table(so@meta.data[[cond_col]][is_healthy]))


# ── COMPARISON 4: Healthy kidney vs bladder (in vivo) ────────────────────────

message("\n\n########## COMPARISON 4: Healthy in vivo — kidney vs bladder ##########")

so_healthy_invivo <- so[, is_healthy & so$Categories %in% c("kidneyMouse", "bladderMouse")]
so_healthy_invivo$plot_group <- so_healthy_invivo$Categories

message(sprintf("  Healthy in vivo subset: %d cells", ncol(so_healthy_invivo)))
print(table(so_healthy_invivo$Categories))

if (all(c("kidneyMouse", "bladderMouse") %in% so_healthy_invivo$Categories)) {
  run_comparison(
    so_healthy_invivo,
    ident1 = "kidneyMouse",
    ident2 = "bladderMouse",
    label  = "Healthy_KidneyInvivo_vs_BladderInvivo"
  )
} else {
  message("  WARNING: one or both healthy in vivo groups are empty — skipping")
}


# ── COMPARISON 5: Healthy organoids — pairwise ───────────────────────────────

message("\n\n########## COMPARISON 5: Healthy organoids — pairwise ##########")

# Organoid samples are inherently healthy/control-derived, but apply the same
# healthy filter in case any disease organoid conditions are present.
so_healthy_org <- so[, is_healthy & so$Categories %in% c("kidneyOrganoid",
                                                          "bladderOrganoid",
                                                          "ureterOrganoid")]
so_healthy_org$plot_group <- so_healthy_org$Categories

present_healthy_org <- intersect(
  c("kidneyOrganoid", "bladderOrganoid", "ureterOrganoid"),
  unique(so_healthy_org$Categories)
)

message(sprintf("  Healthy organoid subset: %d cells", ncol(so_healthy_org)))
print(table(so_healthy_org$Categories))

if (length(present_healthy_org) >= 2) {
  pairs_org <- combn(present_healthy_org, 2, simplify = FALSE)
  for (pair in pairs_org) {
    run_comparison(
      so_healthy_org,
      ident1 = pair[1],
      ident2 = pair[2],
      label  = paste0("Healthy_", pair[1], "_vs_", pair[2])
    )
  }
} else {
  message("  WARNING: fewer than 2 healthy organoid groups found — skipping")
}


################################################################################
# Summary table of all comparison sizes
################################################################################

message("\n\n========== Summary ==========")
summary_df <- data.frame(
  Categories = names(table(so$Categories)),
  n_cells    = as.integer(table(so$Categories))
)
print(summary_df)
write.csv(summary_df, file.path(CMP_DIR, "group_cell_counts.csv"), row.names = FALSE)

message("\nAll comparisons complete.")
message("Outputs in: ", CMP_DIR)






