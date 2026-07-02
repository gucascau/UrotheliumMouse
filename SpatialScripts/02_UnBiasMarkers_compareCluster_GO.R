################################################################################
# 02_UnBiasMarkers_compareCluster_GO.R
#
# DEG detection + comparative functional enrichment for the unbiased
# FindAllMarkers result (UnBiasMarkersSpatialHD, from
# 01_VisiumHD_integrate_harmony_sketch.R, kidney3p + kidney harmony/sketch
# clusters).
#
# Cutoff (chosen after checking per-cluster DEG yield — see note below):
#   p_val_adj < 0.05, avg_log2FC > 0.58 (~1.5-fold), pct.1 > 0.05
#
# NOTE on pct.1 threshold: the conventional pct.1 > 0.1 cutoff was tested
# first and left cluster 5 with only 1 DEG at log2FC > 1 — its high-log2FC
# genes are detected in just 1-2% of cluster-5 bins (sparse-detection
# artifact), too few for meaningful enrichment. Relaxing to pct.1 > 0.05
# brings every cluster to a usable DEG count.
#
# Two complementary enrichment strategies, each for GO Biological Process
# and KEGG:
#   ORA  (over-representation, compareCluster + enrichGO/enrichKEGG) —
#        tests the padj/log2FC/pct.1-filtered DEG list per cluster against
#        the annotated background.
#   GSEA (gene set enrichment, compareCluster + gseGO/gseKEGG) — tests the
#        full per-cluster ranked gene list (all genes FindAllMarkers tested
#        for that cluster, ranked by avg_log2FC, not just the DEG subset),
#        so it also captures coordinated sub-threshold shifts.
#
# Steps:
#   1. Load UnBiasMarkersSpatialHD.rds
#   2. Filter DEGs per cluster (padj < 0.05, avg_log2FC > 0.58, pct.1 > 0.05)
#   3. compareCluster ORA  — GO Biological Process (enrichGO)
#   4. compareCluster ORA  — KEGG (enrichKEGG, Entrez-mapped DEGs)
#   5. compareCluster GSEA — GO Biological Process (gseGO)
#   6. compareCluster GSEA — KEGG (gseKEGG, Entrez-mapped ranked lists)
#
# Output: SpatialScripts/output/GO_enrichment/
#           ├─ DEGs_padj0.05_log2FC0.58_pct0.05.csv
#           ├─ compareCluster_ORA_GO_BP.[rds|csv], _dotplot.pdf
#           ├─ compareCluster_ORA_KEGG.[rds|csv], _dotplot.pdf
#           ├─ compareCluster_GSEA_GO_BP.[rds|csv], _dotplot.pdf
#           └─ compareCluster_GSEA_KEGG.[rds|csv], _dotplot.pdf
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(ggplot2)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
OUT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/SpatialScripts/output"
ENR_DIR <- file.path(OUT_DIR, "GO_enrichment")
dir.create(ENR_DIR, showWarnings = FALSE, recursive = TRUE)

MARKERS_RDS <- file.path(OUT_DIR, "UnBiasMarkersSpatialHD.rds")

# ── Parameters ────────────────────────────────────────────────────────────────
PADJ_THRESH  <- 0.05
LFC_THRESH   <- 0.58
PCT1_THRESH  <- 0.1
KEGG_ORGANISM <- "mmu"
TOP_PATHWAY  <- 10   # categories shown per cluster in the dot-plots

# ── Step 1: Load markers ───────────────────────────────────────────────────────
message("==> Loading: ", MARKERS_RDS)
UnBiasMarkersSpatialHD <- readRDS(MARKERS_RDS)

# ── Step 2: Filter DEGs ────────────────────────────────────────────────────────
message(sprintf("==> Filtering DEGs (padj < %.2f, log2FC > %.2f, pct.1 > %.2f) ...",
                PADJ_THRESH, LFC_THRESH, PCT1_THRESH))

deg_df <- UnBiasMarkersSpatialHD %>%
  filter(p_val_adj < PADJ_THRESH, avg_log2FC > LFC_THRESH, pct.1 > PCT1_THRESH) %>%
  arrange(cluster, desc(avg_log2FC))

message(sprintf("  %d DEGs across %d clusters", nrow(deg_df), n_distinct(deg_df$cluster)))
print(deg_df %>% count(cluster) %>% arrange(as.numeric(as.character(cluster))))

write.csv(deg_df,
  file.path(ENR_DIR, sprintf("DEGs_padj%.2f_log2FC%.2f_pct%.2f.csv",
                              PADJ_THRESH, LFC_THRESH, PCT1_THRESH)),
  row.names = FALSE)

deg_df$cluster <- factor(deg_df$cluster,
  levels = sort(unique(as.numeric(as.character(deg_df$cluster)))))

deg_df %>% as.data.frame()%>% filter(cluster == 9) %>% arrange(desc(avg_log2FC)) %>% head(10) %>%
  dplyr::select(gene, avg_log2FC, pct.1, pct.2, p_val_adj)
# Save + dot-plot helper shared by all four enrichment runs
save_compareCluster <- function(cc_res, tag, title) {
  if (is.null(cc_res) || nrow(cc_res@compareClusterResult) == 0) {
    message(sprintf("  [%s] no enriched terms — skipping save", tag))
    return(invisible(NULL))
  }
  message(sprintf("  [%s] %d result rows across clusters", tag, nrow(cc_res@compareClusterResult)))

  saveRDS(cc_res, file.path(ENR_DIR, sprintf("compareCluster_%s.rds", tag)))
  write.csv(as.data.frame(cc_res),
    file.path(ENR_DIR, sprintf("compareCluster_%s.csv", tag)), row.names = FALSE)

  p <- dotplot(cc_res, showCategory = TOP_PATHWAY) +
    ggtitle(title) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(ENR_DIR, sprintf("compareCluster_%s_dotplot.pdf", tag)),
         plot = p, width = 12, height = 10)
  invisible(cc_res)
}

# ── Step 3: compareCluster ORA — GO Biological Process ────────────────────────
message("==> Running compareCluster ORA (enrichGO, BP) ...")

cc_ora_go <- compareCluster(
  gene          ~ cluster,
  data          = deg_df,
  fun           = "enrichGO",
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 1,
  readable      = TRUE
)
save_compareCluster(cc_ora_go, "ORA_GO_BP",
  "GO Biological Process (ORA) — unbiased cluster DEGs")

# we selected some representative GO terms for the manuscript figure, but the full results are in the output CSV and RDS files.
SelectedPathways<- c("epidermis development","epidermal cell differentiation","morphogenesis of a branching epithelium","kidney epithelium development","keratinocyte differentiation","mesenchymal cell differentiation","epithelial to mesenchymal transition","nephron epithelium development","nephron morphogenesis","Notch signaling pathway","wound healing","T cell receptor signaling pathway","leukocyte migration","macrophage migration","regulation of leukocyte cell-cell adhesion","cell-cell junction organization","muscle cell differentiation","cell-substrate adhesion")
SelectedGODotPlot <- dotplot(cc_ora_go, showCategory = SelectedPathways, label_format=NULL) +
    #ggtitle("GO Biological Process (ORA) — selected pathways") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(ENR_DIR, "compareCluster_selectedpathways_dotplot.pdf"),
         plot = SelectedGODotPlot, width = 8, height = 6)
? dotplot

# ── Step 4: compareCluster ORA — KEGG ──────────────────────────────────────────
message("==> Mapping DEG symbols to Entrez IDs for KEGG ...")

deg_entrez <- bitr(unique(deg_df$gene), fromType = "SYMBOL", toType = "ENTREZID",
                    OrgDb = org.Mm.eg.db) %>%
  distinct(SYMBOL, .keep_all = TRUE)

deg_df_entrez <- deg_df %>%
  inner_join(deg_entrez, by = c("gene" = "SYMBOL")) %>%
  transmute(gene = ENTREZID, cluster)

message(sprintf("  %d / %d DEGs mapped to Entrez", nrow(deg_df_entrez), nrow(deg_df)))

message("==> Running compareCluster ORA (enrichKEGG) ...")
cc_ora_kegg <- compareCluster(
  gene          ~ cluster,
  data          = deg_df_entrez,
  fun           = "enrichKEGG",
  organism      = KEGG_ORGANISM,
  pAdjustMethod = "BH",
  pvalueCutoff  = 1
)
save_compareCluster(cc_ora_kegg, "ORA_KEGG",
  "KEGG pathways (ORA) — unbiased cluster DEGs")

# ── Step 5: compareCluster GSEA — GO Biological Process ───────────────────────
# Ranked lists use ALL genes FindAllMarkers tested per cluster (not just the
# DEG-filtered subset), ranked by avg_log2FC, so GSEA can also pick up
# coordinated sub-threshold shifts that ORA on the strict DEG list would miss.
message("==> Building per-cluster ranked gene lists for GSEA ...")

build_ranked_list <- function(cluster_id, id_col = "gene") {
  sub <- UnBiasMarkersSpatialHD %>%
    filter(cluster == cluster_id) %>%
    filter(!is.na(.data[[id_col]]), !is.na(avg_log2FC)) %>%
    distinct(.data[[id_col]], .keep_all = TRUE) %>%
    arrange(desc(avg_log2FC))
  setNames(sub$avg_log2FC, sub[[id_col]])
}

cluster_ids <- sort(unique(as.numeric(as.character(UnBiasMarkersSpatialHD$cluster))))

gene_lists_symbol <- setNames(
  lapply(cluster_ids, build_ranked_list, id_col = "gene"),
  cluster_ids
)

message("==> Running compareCluster GSEA (gseGO, BP) ...")
cc_gsea_go <- compareCluster(
  geneClusters  = gene_lists_symbol,
  fun           = "gseGO",
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 1,
  seed          = TRUE
)
save_compareCluster(cc_gsea_go, "GSEA_GO_BP",
  "GO Biological Process (GSEA) — unbiased cluster ranked genes")

# ── Step 6: compareCluster GSEA — KEGG ─────────────────────────────────────────
message("==> Mapping full marker table to Entrez IDs for GSEA-KEGG ...")

all_entrez <- bitr(unique(UnBiasMarkersSpatialHD$gene),
                    fromType = "SYMBOL", toType = "ENTREZID",
                    OrgDb = org.Mm.eg.db) %>%
  distinct(SYMBOL, .keep_all = TRUE)

markers_entrez <- UnBiasMarkersSpatialHD %>%
  inner_join(all_entrez, by = c("gene" = "SYMBOL")) %>%
  transmute(cluster, gene = ENTREZID, avg_log2FC)

gene_lists_entrez <- setNames(
  lapply(cluster_ids, function(cl) {
    sub <- markers_entrez %>%
      filter(cluster == cl) %>%
      distinct(gene, .keep_all = TRUE) %>%
      arrange(desc(avg_log2FC))
    setNames(sub$avg_log2FC, sub$gene)
  }),
  cluster_ids
)

message("==> Running compareCluster GSEA (gseKEGG) ...")
cc_gsea_kegg <- compareCluster(
  geneClusters  = gene_lists_entrez,
  fun           = "gseKEGG",
  organism      = KEGG_ORGANISM,
  pAdjustMethod = "BH",
  pvalueCutoff  = 1,
  seed          = TRUE
)
save_compareCluster(cc_gsea_kegg, "GSEA_KEGG",
  "KEGG pathways (GSEA) — unbiased cluster ranked genes")

message("==> Done. Outputs in: ", ENR_DIR)
