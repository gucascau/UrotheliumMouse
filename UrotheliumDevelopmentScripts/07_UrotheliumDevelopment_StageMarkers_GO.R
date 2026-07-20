################################################################################
# 07_UrotheliumDevelopment_StageMarkers_GO.R
#
# Candidate marker discovery per developmental stage (Uro cells only,
# n=715), followed by GO Biological Process enrichment on each stage's
# marker set -- a data-driven complement to the curated Fig2/Fig8 gene
# panels: what actually comes up as differential per stage, unbiased, and
# what biology does that imply.
#
# Markers: Seurat::FindAllMarkers (Wilcoxon), one stage vs. all other
# stages pooled, on the log-normalized "originalexp" assay (no counts layer
# in this object -- see 01/05's notes -- FindMarkers' Wilcoxon test only
# needs the data slot). only.pos = TRUE: "candidate markers for stage X"
# means genes UP in X, not down (down-in-X is just "up in some other
# stage" and shows up there instead). Significant = p_val_adj < 0.05 &
# avg_log2FC > 0.25, same log2FC convention as Seurat's default.
#
# GO enrichment: clusterProfiler::enrichGO (BP ontology) run separately per
# stage's significant marker set, against a detected-gene background (any
# gene with nonzero "data" in >=1 Uro cell) rather than the whole genome --
# standard practice for scRNA marker GO, avoids inflating enrichment from
# genes that were never expressed in this cell type to begin with. Symbol
# -> ENTREZID via org.Mm.eg.db (some markers won't map; dropped with a
# warning count, not silently). Per-stage enrichGO results are redundancy-
# trimmed (clusterProfiler::simplify, semantic similarity cutoff 0.7) then
# combined with merge_result() for one side-by-side comparison dotplot.
#
# Input:  UrotheliumDevelopmentScripts/output/UrotheliumOnly_full_annotated.rds
#         (written by 05_UrotheliumDevelopment_ModuleScores.R)
# Output: UrotheliumDevelopment_StageMarkers_all.csv       (full marker table)
#         UrotheliumDevelopment_StageMarkers_top20.csv     (top 20 per stage)
#         UrotheliumDevelopment_GO_BP_<stage>.csv           (one per stage)
#         Fig11_UrotheliumDevelopment_StageMarkers_DotPlot.pdf
#         Fig12_UrotheliumDevelopment_GO_BP_Comparison.pdf
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(viridisLite)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(pheatmap)
  library(RColorBrewer)
  # dplyr loaded LAST: clusterProfiler/org.Mm.eg.db/AnnotationDbi define S4
  # generics (select, rename, filter, mutate, arrange, slice, group_by, ...)
  # that mask their dplyr namesakes if dplyr loads first, breaking every
  # %>% pipe below with cryptic "object not found" errors. Loading dplyr
  # last puts it at the top of the search path so its verbs win.
  library(dplyr)
})

SCRIPT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/UrotheliumDevelopmentScripts"
OUT_DIR    <- file.path(SCRIPT_DIR, "output")

STAGE_ORDER  <- c("E16.5", "P0", "W3", "W12", "W52", "W92")
STAGE_COLORS <- setNames(viridisLite::viridis(length(STAGE_ORDER)), STAGE_ORDER)

N_TOP_MARKERS_PLOT <- 6      # markers per stage shown in Fig11
N_TOP_GO_PLOT       <- 6      # GO BP terms per stage shown in Fig12
PADJ_MARKER_CUTOFF  <- 0.05
LOGFC_MARKER_CUTOFF <- 0.25

# ── Load ─────────────────────────────────────────────────────────────────────
message("==> Loading Uro-only annotated object ...")
in_rds <- file.path(OUT_DIR, "UrotheliumOnly_full_annotated.rds")
if (!file.exists(in_rds)) {
  stop("Missing ", in_rds, " -- run 05_UrotheliumDevelopment_ModuleScores.R first.")
}
uro <- readRDS(in_rds)
uro$Age <- factor(uro$Age, levels = STAGE_ORDER)
DefaultAssay(uro) <- "originalexp"
Idents(uro) <- uro$Age
message("  Uro cells per stage:")
print(table(uro$Age))

fm <- uro[["originalexp"]]@meta.features

################################################################################
# Step 1: candidate markers per stage (FindAllMarkers)
################################################################################
message("\n==> Finding candidate markers per stage (Wilcoxon, one-vs-rest) ...")

markers <- FindAllMarkers(uro, assay = "originalexp", slot = "data",
                           only.pos = TRUE,
                           logfc.threshold = 0,
                           test.use = "MAST")

markers$gene_symbol <- fm$gene_symbols[match(markers$gene, rownames(fm))]
markers <- markers %>%
  rename(ensembl_id = gene, stage = cluster) %>%
  mutate(stage = factor(stage, levels = STAGE_ORDER)) %>%
  arrange(stage, p_val_adj, desc(avg_log2FC)) %>%
  # dplyr::select explicit: org.Mm.eg.db/AnnotationDbi (loaded above) define
  # an S4 "select" generic that masks dplyr::select on the search path,
  # breaking this pipe with a cryptic "object 'gene' not found" error.
  dplyr::select(stage, gene_symbol, ensembl_id, avg_log2FC, pct.1, pct.2, p_val, p_val_adj)

message("  Marker counts per stage (p_val_adj < ", PADJ_MARKER_CUTOFF, "):")
markers_sig <- markers %>% filter(p_val_adj < PADJ_MARKER_CUTOFF & avg_log2FC > 0.25, pct.1 >0.1)
print(table(markers_sig$stage))
saveRDS(markers, file.path(OUT_DIR, "UrotheliumDevelopment_StageMarkers_all.rds"))
write.csv(markers, file.path(OUT_DIR, "UrotheliumDevelopment_StageMarkers_all.csv"),
          row.names = FALSE)
message("  Saved: UrotheliumDevelopment_StageMarkers_all.csv")

top20 <- markers_sig %>% group_by(stage) %>% slice_head(n = 20) %>% ungroup()
write.csv(top20, file.path(OUT_DIR, "UrotheliumDevelopment_StageMarkers_top20.csv"),
          row.names = FALSE)
message("  Saved: UrotheliumDevelopment_StageMarkers_top20.csv")
markers_sig %>% head()
################################################################################
# Figure 11: top marker dot plot by stage
################################################################################
message("\n==> Building Figure 11 (stage marker dot plot) ...")

top_plot <- markers_sig %>% filter(!grepl("mt", gene_symbol)) %>%
  group_by(stage) %>%
  slice_head(n = N_TOP_MARKERS_PLOT) %>%
  ungroup()

# Symbol-keyed "Markers" assay for plotting (rownames(uro) are Ensembl IDs
# in the full object) -- same pattern as 01's Figure 2.
plot_ens <- unique(top_plot$ensembl_id)
mat <- GetAssayData(uro, assay = "originalexp", layer = "data")[plot_ens, , drop = FALSE]
rownames(mat) <- fm$gene_symbols[match(plot_ens, rownames(fm))]
uro[["Markers"]] <- CreateAssayObject(data = mat)
DefaultAssay(uro) <- "Markers"

# Feature order: grouped by stage (left-to-right = E16.5 -> W92), duplicates
# (a gene appearing as top marker in >1 stage) dropped after first occurrence
# so each gene is plotted once.
feature_order <- top_plot %>% distinct(gene_symbol, .keep_all = TRUE) %>% pull(gene_symbol)

fig11 <- DotPlot(uro, features = feature_order, group.by = "Age", assay = "Markers",
                  cols = c("lightgrey", "red"), dot.scale = 7) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.title  = element_text(face = "bold", size = 14)) +
  guides(color = guide_colorbar(title = "AveExp"), size = guide_legend(title = "PerExp")) +
  labs(x = NULL, y = "Developmental stage",
       title = sprintf("Figure 11. Top %d candidate markers per stage", N_TOP_MARKERS_PLOT),
       subtitle = "One-vs-rest MAST, p_adj < 0.05, ranked by avg log2FC")

ggsave(file.path(OUT_DIR, "Fig11_UrotheliumDevelopment_StageMarkers_DotPlot.pdf"),
       fig11, width = max(9, 0.35 * length(feature_order)), height = 4.5)
message("  Saved: Fig11_UrotheliumDevelopment_StageMarkers_DotPlot.pdf")

DefaultAssay(uro) <- "originalexp"

################################################################################
# Step 2: GO Biological Process enrichment per stage
################################################################################
message("\n==> Mapping marker symbols to ENTREZID (org.Mm.eg.db) ...")

# Detected-gene background: any gene with nonzero expression in >=1 Uro
# cell, mapped to ENTREZID -- the enrichment universe, not the whole genome.
detected_ens <- rownames(uro)[Matrix::rowSums(GetAssayData(uro, assay = "originalexp",
                                                             layer = "data") > 0) > 0]
detected_symbols <- unique(na.omit(fm$gene_symbols[match(detected_ens, rownames(fm))]))

bitr_map <- suppressMessages(
  bitr(detected_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
)
message(sprintf("  %d / %d detected gene symbols mapped to ENTREZID (background universe)",
                 nrow(bitr_map), length(detected_symbols)))
universe_entrez <- unique(bitr_map$ENTREZID)

markers_sig <- markers_sig %>%
  left_join(bitr_map, by = c("gene_symbol" = "SYMBOL"))
n_unmapped <- sum(is.na(markers_sig$ENTREZID))
if (n_unmapped > 0) {
  message(sprintf("  %d / %d significant markers had no ENTREZID mapping (dropped from GO input)",
                   n_unmapped, nrow(markers_sig)))
}

message("\n==> Running GO BP enrichment across stages (clusterProfiler::compareCluster) ...")

# Per-stage gene lists (>=5 mappable ENTREZID genes, same floor as before);
# split(..., drop = TRUE) removes stages with zero genes left after the
# ENTREZID join so compareCluster never sees an empty cluster.
stage_gene_counts <- table(markers_sig$stage[!is.na(markers_sig$ENTREZID)])
for (st in STAGE_ORDER) {
  message(sprintf("  %-6s: %d marker genes -> %d with ENTREZID", st,
                   sum(markers_sig$stage == st), stage_gene_counts[st]))
}
keep_stages <- names(stage_gene_counts)[stage_gene_counts >= 5]
if (length(keep_stages) < length(STAGE_ORDER)) {
  message("  Excluding stage(s) with <5 mappable marker genes: ",
          paste(setdiff(STAGE_ORDER, keep_stages), collapse = ", "))
}

# we only used the upregulated genes in each stage
gene_clusters <- markers_sig %>% filter(!grepl("mt-", gene_symbol)) %>%
  filter(!is.na(ENTREZID), stage %in% keep_stages) %>%
  mutate(stage = factor(stage, levels = keep_stages))
gene_clusters %>% tail()

gene_list <- split(gene_clusters$ENTREZID, gene_clusters$stage, drop = TRUE)

################################################################################
# Figure 12: GO BP comparison across stages
################################################################################

# Run the GO enrichment 
cc <- compareCluster(geneCluster = gene_list, fun = "enrichGO",
                        OrgDb = org.Mm.eg.db, ont = "BP", universe = universe_entrez,
                        pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                        readable = TRUE)
# save the RDS file
saveRDS(cc, file = file.path(OUT_DIR, "UrothliumDevelopment_Stage_GOenrichment.rds"))

write.csv(cc, file.path(OUT_DIR, "UrothliumDevelopment_Stage_GOenrichment.csv "))
# for the selected pathways
SelectedGOTerms<- c("stem cell differentiation","embryonic organ development","post-embryonic development","mesenchyme development","morphogenesis of an epithelium","epithelial cell migration","epithelial cell proliferation","metanephros development","mesonephric tubule development","ureteric bud development","epithelium migration","Wnt signaling pathway","non-canonical Wnt signaling pathway","insulin receptor signaling pathway","insulin-like growth factor receptor signaling pathway","angiogenesis","synapse organization","focal adhesion assembly","collagen-activated signaling pathway" )
fig12 <- dotplot(cc, showCategory = SelectedGOTerms, font.size = 6.5) +RotatedAxis()
ggsave(file.path(OUT_DIR, "Fig12_UrotheliumDevelopment_GO_BP_Comparison.pdf"),
           fig12, width = 5, height = 6.5)

################################################################################
# Figure 13: stem cell differentiation / post-embryonic development genes,
# heatmap across developmental stages
################################################################################
message("\n==> Building Figure 13 (stem cell diff. / post-embryonic dev. gene heatmap) ...")

HEATMAP_GO_TERMS <- c("stem cell differentiation", "post-embryonic development", "mesenchyme development")
go_hits <- as.data.frame(cc) %>% filter(Description %in% HEATMAP_GO_TERMS)

if (nrow(go_hits) == 0) {
  message("  Neither GO term found in the compareCluster results -- skipping Figure 13.")
} else {
  print(go_hits %>% dplyr::select(Cluster, Description, p.adjust, geneID))

  # Gene -> GO term membership (a gene hit by both terms gets a combined
  # label), used for the heatmap's row annotation.
  gene_term_map <- do.call(rbind, lapply(seq_len(nrow(go_hits)), function(i) {
    data.frame(gene_symbol = strsplit(go_hits$geneID[i], "/")[[1]],
               GOterm = go_hits$Description[i], stringsAsFactors = FALSE)
  })) %>%
    group_by(gene_symbol) %>%
    summarise(GOterm = paste(unique(GOterm), collapse = " & "), .groups = "drop")

  heatmap_genes <- gene_term_map$gene_symbol
  message(sprintf("  %d unique gene(s) from %d GO term hit(s): %s",
                   length(heatmap_genes), nrow(go_hits), paste(heatmap_genes, collapse = ", ")))

  heatmap_ens <- rownames(fm)[match(heatmap_genes, fm$gene_symbols)]

  # Per-stage mean of the log-normalized "data" values (same convention as
  # 05's module scores -- no un-logging/re-logging), one column per
  # developmental stage (a "condition"), not one column per cell.
  mat_hm <- GetAssayData(uro, assay = "originalexp", layer = "data")[heatmap_ens, , drop = FALSE]
  rownames(mat_hm) <- heatmap_genes

  avg_exp <- sapply(STAGE_ORDER, function(st) {
    Matrix::rowMeans(mat_hm[, uro$Age == st, drop = FALSE])
  })
  rownames(avg_exp) <- heatmap_genes

  annotation_row <- data.frame(GOterm = gene_term_map$GOterm, row.names = gene_term_map$gene_symbol)
  row_annot_colors <- setNames(
    c("#D55E00", "#0072B2", "#009E73")[seq_along(unique(annotation_row$GOterm))],
    unique(annotation_row$GOterm)
  )
  annotation_col <- data.frame(Stage = STAGE_ORDER, row.names = STAGE_ORDER)

  # Diverging blue-white-red scale (classic publication convention for
  # row-scaled/z-scored expression, not a sequential palette like viridis --
  # white sits at z = 0 so "above/below the gene's own mean across stages"
  # reads at a glance). Breaks fixed to +/-2.5 SD (not data-driven) so color
  # is comparable across re-runs and a couple of outlier genes can't
  # compress the rest of the scale.
  hm_palette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
  hm_breaks  <- seq(-2.5, 2.5, length.out = length(hm_palette) + 1)

  pdf(file.path(OUT_DIR, "Fig13_UrotheliumDevelopment_StemCell_PostEmbryonic_Heatmap.pdf"),
      width = 4, height =6)
  pheatmap(avg_exp, scale = "row", cluster_cols = FALSE, cluster_rows = TRUE,
           # Reverses the row dendrogram's leaf order (top<->bottom) without
           # changing the clustering itself -- pheatmap has no plain
           # "reverse" argument, so this is the standard way via hclust.
           clustering_callback = function(hc, ...) as.hclust(rev(as.dendrogram(hc))),
           color = hm_palette, breaks = hm_breaks, border_color = NA,
           #annotation_row = annotation_row, annotation_col = annotation_col,
           #annotation_colors = list(GOterm = row_annot_colors, Stage = STAGE_COLORS),
           main = "Stem cell differentiation / post-embryonic development genes")
  dev.off()
  message("  Saved: Fig13_UrotheliumDevelopment_StemCell_PostEmbryonic_Heatmap.pdf")
}

message("\n==> Done.")
