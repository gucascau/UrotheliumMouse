#!/usr/bin/env Rscript
################################################################################
# 02_UrotheliumDevelopment_Spp1_CellChat.R
#
# Cell-cell communication across urothelium development, focused on the Spp1
# (osteopontin) signaling pathway: is Urothelium a source and/or target of
# Spp1 signal, and how does that change from E16.5 through W92?
#
# Motivation: UrotheliumDevelopmentScripts/05_UrotheliumDevelopment_
# ModuleScores.R included Spp1 in an exploratory "Stress / remodeling"
# module scored ON Uro cells (i.e. Spp1 expression within Uro itself); the
# separate CellCommunicationScripts/01_Urothelium_CellChat_Communication.R
# run (different dataset -- postnatal healthy-kidney atlas, not this
# developmental series) also flagged Spp1 as a top Urothelium-sender ligand.
# This script asks the actual cell-cell communication question -- who is
# Urothelium signaling to/from via Spp1 -- on the developmental atlas
# specifically, resolved by stage.
#
# Design (per user's explicit choices):
#   - One CellChat run PER developmental stage (E16.5/P0/W3/W12/W52/W92),
#     not pooled -- lets Spp1 signaling strength/targets be tracked across
#     development, same spirit as Fig8's module-score-by-Age trend.
#   - Full celltype_final (26 types) as the group.by grouping (same
#     granularity as 01_Urothelium_CellChat_Communication.R), not the
#     collapsed 4-class CellClass from 01_UrotheliumDevelopment_Figures.R --
#     resolves which specific cell type(s) send/receive Spp1, at the cost of
#     sparse counts for rare types at some stages (e.g. NP/IM vanish after
#     W3 -- expected developmental biology, not a bug; celltype x stage
#     counts checked interactively before writing this script).
#   - Full CellChatDB.mouse run per stage (not restricted to Spp1 upfront)
#     so per-stage global outputs are available too, then every downstream
#     step subsets to the SPP1 pathway specifically. Matches
#     01_Urothelium_CellChat_Communication.R's "full DB, then focus"
#     pattern. CellChatDB.mouse's SPP1 pathway = 9 L-R pairs: SPP1_CD44 and
#     8 SPP1_ITGA*_ITGB* integrin heterodimers (checked interactively).
#
# Input data: the converted Chen2025 RDS carries only a log-normalized
# "data" layer (no counts, no scale.data) -- see 01_UrotheliumDevelopment_
# Figures.R's header notes. CellChat's createCellChat()/computeCommunProb()
# operate on normalized expression, not raw counts (same assumption
# 01_Urothelium_CellChat_Communication.R already makes for its input
# object), so this log-normalized "data" is used directly -- the raw-count
# h5ad export (VisiumLowDevelopmentScripts/output/Chen2025_rawcounts, built
# for RCTD's Poisson model) is not needed here.
#
# rownames(obj) are Ensembl IDs; CellChatDB.mouse keys on gene symbols. All
# 31,671 gene_symbols in meta.features are unique and non-missing (checked
# interactively), so the whole assay is safely renamed to symbols once,
# per-stage (after subsetting, to keep peak memory down), rather than
# building a small marker-panel assay as 01/05 do.
#
# Input:  RawMouseSingleCellDatasets/MultiOmicSpatialMouseKidney_SexSpecific_
#         lifespan_Chen2025NatGenet_zellkonvertedConverted.rds
# Output: output/UroDev_Spp1/CellChat_UroDev_<stage>.rds   (per-stage CellChat objects)
#         output/UroDev_Spp1/SPP1_LR_interactions_<stage>.csv
#         output/UroDev_Spp1/SPP1_LR_interactions_AllStages.csv (combined, long format)
#         output/UroDev_Spp1/Spp1_Uro_AsSender_Bubble_<stage>.pdf
#         output/UroDev_Spp1/Spp1_Uro_AsReceiver_Bubble_<stage>.pdf
#         output/UroDev_Spp1/Fig_Spp1_UroDev_SenderTrend_Heatmap.pdf
#         output/UroDev_Spp1/Fig_Spp1_UroDev_ReceiverTrend_Heatmap.pdf
#         output/UroDev_Spp1/DotPlot_Uro_Spp1_Col4a345_ByStage.pdf
#         output/UroDev_Spp1/VlnPlot_Uro_Spp1_Col4a345_ByStage.pdf
#         output/UroDev_Spp1/DotPlot_Receptors_OtherCellTypes_ByAge.pdf
################################################################################

suppressMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(viridisLite)
  library(CellChat)
  library(future)
})

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 2e11)
options(future.rng.onMisuse = "ignore")
set.seed(200000)

n_workers <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))
plan("multicore", workers = n_workers)

# ── Paths ─────────────────────────────────────────────────────────────────────
RAW_BASE   <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/RawMouseSingleCellDatasets"
SCRIPT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/CellCommunicationScripts"
OUT_DIR    <- file.path(SCRIPT_DIR, "output", "UroDev_Spp1")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OBJECT_RDS <- file.path(RAW_BASE,
  "MultiOmicSpatialMouseKidney_SexSpecific_lifespan_Chen2025NatGenet_zellkonvertedConverted.rds")

STAGE_ORDER  <- c("E16.5", "P0", "W3", "W12", "W52", "W92")
STAGE_COLORS <- setNames(viridisLite::viridis(length(STAGE_ORDER)), STAGE_ORDER)
UroLabel     <- "Uro"
SPP1_PATHWAY <- "SPP1"
MIN_CELLS_PER_GROUP <- 10   # groups below this are dropped per-stage before CellChat

# ── Load full developmental atlas ───────────────────────────────────────────
message("==> Loading Chen2025 developmental atlas object ...")
obj <- readRDS(OBJECT_RDS)
message(sprintf("  %d cells x %d genes", ncol(obj), nrow(obj)))
DefaultAssay(obj) <- "originalexp"
obj$Age <- factor(obj$Age, levels = STAGE_ORDER)
obj$celltype_final <- droplevels(factor(obj$celltype_final))
fm <- obj[["originalexp"]]@meta.features
stopifnot(!any(is.na(fm$gene_symbols)), !any(duplicated(fm$gene_symbols)))

CellChatDB <- CellChatDB.mouse

################################################################################
# Uro-only expression: Spp1 (ligand) vs Col4a3/Col4a4/Col4a5 (BM collagen IV
# chains), by stage
#
# Purely descriptive -- no CellChat network math, no per-stage cell-type
# filtering -- restricted to Uro cells across all 6 stages at once. Spp1 is
# the ligand this whole script tracks via CellChat; Col4a3/Col4a4/Col4a5
# form the alpha3-4-5 collagen IV heterotrimer that is the basement-membrane
# scaffold Spp1/integrin signaling acts within. Plotting both in Uro cells
# specifically distinguishes two different questions the CellChat bubble
# plots above can't: is Uro itself ramping up Spp1 ligand production over
# development (ligand induction), vs. does Uro also express the BM collagen
# genes that would indicate it's building/maintaining that niche locally
# rather than just receiving it from stroma (receptor/niche availability).
################################################################################
message("\n==> Uro-only Spp1/Col4a3/Col4a4/Col4a5 expression across stages ...")

URO_GENES <- c("Spp1", "Col4a3", "Col4a4", "Col4a5", "App")
uro_gene_ids <- setNames(rownames(fm)[match(URO_GENES, fm$gene_symbols)], URO_GENES)
missing_uro_genes <- URO_GENES[is.na(uro_gene_ids)]

uro_obj <- subset(obj, subset = celltype_final == UroLabel)
# uro_obj <- subset(uro_obj, features = uro_gene_ids)
message(sprintf("  %d Uro cells across stages: %s", ncol(uro_obj),
                paste(table(uro_obj$Age), collapse = ", ")))

# rownames(uro_obj) is intentionally left as Ensembl IDs -- Seurat silently
# no-ops renaming features on this object's assay ("Renaming features in
# v3/v4 assays is not supported", confirmed interactively), so relabeling to
# gene symbols happens on the plot axes below instead of on the object.
symbol_by_id <- setNames(names(uro_gene_ids), uro_gene_ids)


DotPlot_Uro_Spp1Col4 <- DotPlot(
    uro_obj,
    features  = rev(unname(uro_gene_ids)),
    group.by  = "Age",
    cols      = c("lightgrey", "red"),
    dot.scale = 8
  ) +
    scale_x_discrete(labels = symbol_by_id[unname(uro_gene_ids)]) +  RotatedAxis()  + coord_flip()+ 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "left") +
    guides(color = guide_colorbar(title = "AveExp"), size = guide_legend(title = "PerExp"))

  ggsave(file.path(OUT_DIR, "DotPlot_Uro_Spp1_Col4a345_ByStage.pdf"),
         DotPlot_Uro_Spp1Col4, width = 4, height = 4.5)
  message("  Saved: DotPlot_Uro_Spp1_Col4a345_ByStage.pdf")

  vln_plots <- lapply(unname(uro_gene_ids), function(gid) {
    VlnPlot(uro_obj, features = gid, group.by = "Age", cols = STAGE_COLORS, pt.size = 0) +
      ggtitle(symbol_by_id[[gid]]) +
      theme(legend.position = "none", axis.title.x = element_blank())
  })
  VlnPlot_Uro_Spp1Col4 <- patchwork::wrap_plots(vln_plots, ncol = 2)
  ggsave(file.path(OUT_DIR, "VlnPlot_Uro_Spp1_Col4a345_ByStage.pdf"),
         VlnPlot_Uro_Spp1Col4, width = 9, height = 7)
  message("  Saved: VlnPlot_Uro_Spp1_Col4a345_ByStage.pdf")

# Receptor dotplot: Spp1-relevant receptors (Cd44 + Itga*/Itgb* integrin
# heterodimer partners forming the CellChatDB SPP1 pathway, plus Sdc4/
# Adgrg6/Cd74/Tnfrsf21) in the non-Uro cell types that could plausibly be
# on the receiving end of Uro-derived Spp1, split by stage. Complements the
# Uro-only Spp1/Col4a3-5 DotPlot above: that one asks whether Uro is
# inducing the ligand, this one asks whether the putative receiving niche
# cell types carry the receptor at each stage.
#
# Built via per-stage DotPlot() + facet_wrap(), NOT DotPlot(split.by="Age")
# -- confirmed interactively that split.by mode has no real AveExp color
# scale to show: Seurat blends each split group's base color (cols) with
# expression via lightness through scale_color_identity(), so there's no
# continuous quantity for guide_colorbar() to draw (hence its "needs
# continuous scales" warning and the missing legend), and worse, each age
# panel would be normalized to its own hue as "max", making color intensity
# incomparable across ages. Faceting on Age with one shared continuous
# AveExp gradient fixes both: a real legend, and cross-age comparability.
SenderReceptors <- c("Cd44","Itgav","Itgb1","Itgb3","Itgb4","Itgb5","Itgb6","Itgb8","Sdc4","Adgrg6","Cd74","Tnfrsf21")

# for the selected cell types -- NP/UBP/IM (nephron progenitor, ureteric bud
# progenitor, intermediate mesoderm) are real populations only at E16.5/P0;
# they're developmentally expected to vanish after that (same "NP/IM vanish
# after W3" point made in this script's header), so they're kept in
# OtherCellTypes but their data is dropped at W3/W12/W52/W92 below, and
# facet_wrap(scales="free_x") (not the default shared x-scale) is what
# actually removes their axis tick from those later-stage panels -- with a
# shared discrete x-scale, dropping the data rows alone still leaves an
# empty tick/label per panel, since all panels share one axis built from
# the full factor level set (confirmed interactively).
obj@meta.data$celltype_final %>% table()
OtherCellTypes <- c("NP", "UBP", "IM", "Endo", "Fib", "Macro", "T", "B/Plasma")
Idents(obj) <- "celltype_final"
OtherObj <- subset(obj, subset = celltype_final %in% OtherCellTypes)
Idents(OtherObj) <- factor(Idents(OtherObj), levels = OtherCellTypes)

receptor_gene_ids <- setNames(rownames(fm)[match(SenderReceptors, fm$gene_symbols)], SenderReceptors)
missing_receptor_genes <- SenderReceptors[is.na(receptor_gene_ids)]
if (length(missing_receptor_genes) > 0) {
  warning("Receptor gene(s) not found in object, dropped: ", paste(missing_receptor_genes, collapse = ", "))
}
receptor_gene_ids <- receptor_gene_ids[!is.na(receptor_gene_ids)]

symbol_by_id_receptor <- setNames(names(receptor_gene_ids), receptor_gene_ids)

receptor_dotplot_data <- bind_rows(lapply(STAGE_ORDER, function(st) {
  stage_cells <- WhichCells(OtherObj, expression = Age == st)
  if (length(stage_cells) == 0) return(NULL)
  stage_obj <- subset(OtherObj, cells = stage_cells)
  d <- DotPlot(stage_obj, features = unname(receptor_gene_ids),
               group.by = "celltype_final")$data
  d$Age <- st
  d
}))
receptor_dotplot_data$features.plot <- factor(receptor_dotplot_data$features.plot,
                                               levels = rev(unname(receptor_gene_ids)))
receptor_dotplot_data$id  <- factor(receptor_dotplot_data$id, levels = OtherCellTypes)
receptor_dotplot_data$Age <- factor(receptor_dotplot_data$Age, levels = STAGE_ORDER)

progenitor_types_early_only <- c("NP", "UBP", "IM")
progenitor_drop_stages <- c("W3", "W12", "W52", "W92")
receptor_dotplot_data <- receptor_dotplot_data %>%
  filter(!(id %in% progenitor_types_early_only & Age %in% progenitor_drop_stages))

DotPlot_receptors <- ggplot(receptor_dotplot_data,
    aes(x = id, y = features.plot, color = avg.exp.scaled, size = pct.exp)) +
  geom_point() +
  scale_color_gradient(low = "lightgrey", high = "red", name = "AveExp") +
  scale_size(range = c(0, 8), name = "PerExp") +
  scale_y_discrete(labels = symbol_by_id_receptor[levels(receptor_dotplot_data$features.plot)]) +
  facet_wrap(~ Age, nrow = 1, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = NULL)

ggsave(file.path(OUT_DIR, "DotPlot_Receptors_OtherCellTypes_ByAge.pdf"),
       DotPlot_receptors,
       width = 10.5,
       height = 4.5)
message("  Saved: DotPlot_Receptors_OtherCellTypes_ByAge.pdf")


################################################################################
# Per-stage CellChat runs
################################################################################
spp1_tables <- list()

for (st in STAGE_ORDER) {
  message(sprintf("\n================ Stage: %s ================", st))

  stage_obj <- subset(obj, subset = Age == st)
  stage_obj$celltype_final <- droplevels(stage_obj$celltype_final)

  ct_counts <- table(stage_obj$celltype_final)
  keep_types <- names(ct_counts)[ct_counts >= MIN_CELLS_PER_GROUP]
  dropped_types <- setdiff(names(ct_counts), keep_types)
  if (length(dropped_types) > 0) {
    message(sprintf("  Dropping %d cell type(s) with <%d cells at %s: %s",
                     length(dropped_types), MIN_CELLS_PER_GROUP, st,
                     paste(dropped_types, collapse = ", ")))
  }
  if (!(UroLabel %in% keep_types)) {
    message(sprintf("  Uro has <%d cells at %s -- skipping stage.", MIN_CELLS_PER_GROUP, st))
    next
  }
  stage_obj <- subset(stage_obj, subset = celltype_final %in% keep_types)
  stage_obj$celltype_final <- droplevels(stage_obj$celltype_final)
  message(sprintf("  %d cells, %d cell types (incl. Uro n=%d)",
                   ncol(stage_obj), nlevels(stage_obj$celltype_final),
                   sum(stage_obj$celltype_final == UroLabel)))

  # Symbol-keyed expression matrix for this stage only (memory-bounded by
  # per-stage cell count, not the full 203k-cell object).
  mat <- GetAssayData(stage_obj, assay = "originalexp", layer = "data")
  rownames(mat) <- fm$gene_symbols[match(rownames(mat), rownames(fm))]

  meta <- stage_obj@meta.data
  meta$celltype_final <- droplevels(meta$celltype_final)

  # ── CellChat object ────────────────────────────────────────────────────────
  cellchat <- createCellChat(object = mat, meta = meta, group.by = "celltype_final")
  cellchat@DB <- CellChatDB

  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 3)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- tryCatch(netAnalysis_computeCentrality(cellchat, slot.name = "netP"),
                        error = function(e) {
                          message("  netAnalysis_computeCentrality failed (", conditionMessage(e), "); continuing without it.")
                          cellchat
                        })

  saveRDS(cellchat, file.path(OUT_DIR, sprintf("CellChat_UroDev_%s.rds", st)))
  message(sprintf("  Saved: CellChat_UroDev_%s.rds", st))

  detected_pathways <- cellchat@netP$pathways
  if (!(SPP1_PATHWAY %in% detected_pathways)) {
    message(sprintf("  SPP1 pathway not detected at %s (no significant L-R pairs passed filtering).", st))
    next
  }

  df_spp1 <- subsetCommunication(cellchat, signaling = SPP1_PATHWAY) %>%
    mutate(stage = st)
  write.csv(df_spp1, file.path(OUT_DIR, sprintf("SPP1_LR_interactions_%s.csv", st)),
            row.names = FALSE)
  message(sprintf("  Saved: SPP1_LR_interactions_%s.csv (%d L-R rows)", st, nrow(df_spp1)))
  spp1_tables[[st]] <- df_spp1

  other_types <- setdiff(levels(cellchat@idents), UroLabel)

  # Uro as Spp1 sender
  spp1_sender <- df_spp1 %>% filter(source == UroLabel)
  if (nrow(spp1_sender) > 0) {
    p_sender <- tryCatch(
      netVisual_bubble(cellchat, sources.use = UroLabel, targets.use = other_types,
                        signaling = SPP1_PATHWAY, remove.isolate = TRUE),
      error = function(e) NULL)
    if (!is.null(p_sender)) {
      ggsave(file.path(OUT_DIR, sprintf("Spp1_Uro_AsSender_Bubble_%s.pdf", st)),
             p_sender, width = 7, height = max(4, length(unique(spp1_sender$interaction_name_2)) / 2))
    }
  } else {
    message(sprintf("  No Uro-as-sender SPP1 interactions passed filtering at %s.", st))
  }

  # Uro as Spp1 receiver
  spp1_receiver <- df_spp1 %>% filter(target == UroLabel)
  if (nrow(spp1_receiver) > 0) {
    p_receiver <- tryCatch(
      netVisual_bubble(cellchat, sources.use = other_types, targets.use = UroLabel,
                        signaling = SPP1_PATHWAY, remove.isolate = TRUE),
      error = function(e) NULL)
    if (!is.null(p_receiver)) {
      ggsave(file.path(OUT_DIR, sprintf("Spp1_Uro_AsReceiver_Bubble_%s.pdf", st)),
             p_receiver, width = 7, height = max(4, length(unique(spp1_receiver$interaction_name_2)) / 2))
    }
  } else {
    message(sprintf("  No Uro-as-receiver SPP1 interactions passed filtering at %s.", st))
  }
}

################################################################################
# Cross-stage synthesis
################################################################################
message("\n==> Combining SPP1 interactions across stages ...")

if (length(spp1_tables) == 0) {
  message("  No stage produced a detectable SPP1 pathway -- nothing to synthesize.")
  quit(save = "no", status = 0)
}

spp1_all <- bind_rows(spp1_tables) %>%
  mutate(stage = factor(stage, levels = STAGE_ORDER))
write.csv(spp1_all, file.path(OUT_DIR, "SPP1_LR_interactions_AllStages.csv"), row.names = FALSE)
message("  Saved: SPP1_LR_interactions_AllStages.csv")

# Uro-as-sender: total Spp1 communication probability to each other cell
# type, by stage (0 where a cell type wasn't present/didn't pass filtering
# at that stage -- absence is itself part of the developmental story).
sender_df <- spp1_all %>%
  filter(source == UroLabel) %>%
  group_by(stage, target) %>%
  summarise(total_prob = sum(prob), .groups = "drop")

if (nrow(sender_df) > 0) {
  fig_sender <- ggplot(sender_df, aes(x = stage, y = target, fill = total_prob)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(name = "Total\nSpp1 prob", option = "viridis", na.value = "grey95") +
    labs(x = "Developmental stage", y = "Target cell type (celltype_final)",
         title = "Spp1 signaling: Urothelium -> other cell types across development",
         subtitle = "Summed CellChat communication probability, SPP1 pathway, per stage") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 13),
          axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(OUT_DIR, "Fig_Spp1_UroDev_SenderTrend_Heatmap.pdf"),
         fig_sender, width = 8, height = max(4, 0.3 * length(unique(sender_df$target))))
  message("  Saved: Fig_Spp1_UroDev_SenderTrend_Heatmap.pdf")
} else {
  message("  No Uro-as-sender SPP1 interactions across any stage -- skipping sender trend figure.")
}

# Uro-as-receiver: total Spp1 communication probability from each other
# cell type into Uro, by stage.
receiver_df <- spp1_all %>%
  filter(target == UroLabel) %>%
  group_by(stage, source) %>%
  summarise(total_prob = sum(prob), .groups = "drop")

if (nrow(receiver_df) > 0) {
  fig_receiver <- ggplot(receiver_df, aes(x = stage, y = source, fill = total_prob)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(name = "Total\nSpp1 prob", option = "viridis", na.value = "grey95") +
    labs(x = "Developmental stage", y = "Source cell type (celltype_final)",
         title = "Spp1 signaling: other cell types -> Urothelium across development",
         subtitle = "Summed CellChat communication probability, SPP1 pathway, per stage") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 13),
          axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(OUT_DIR, "Fig_Spp1_UroDev_ReceiverTrend_Heatmap.pdf"),
         fig_receiver, width = 8, height = max(4, 0.3 * length(unique(receiver_df$source))))
  message("  Saved: Fig_Spp1_UroDev_ReceiverTrend_Heatmap.pdf")
} else {
  message("  No other-cell-type-to-Uro SPP1 interactions across any stage -- skipping receiver trend figure.")
}

message("\n==> Done.")
