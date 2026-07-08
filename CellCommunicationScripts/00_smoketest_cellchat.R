#!/usr/bin/env Rscript
# Smoke test for 01_Urothelium_CellChat_Communication.R
#
# Runs the full CellChat pipeline (createCellChat -> subsetData ->
# identifyOverExpressedGenes/Interactions -> computeCommunProb ->
# filterCommunication -> computeCommunProbPathway -> aggregateNet ->
# netAnalysis_computeCentrality -> Urothelium-focused circle/bubble/
# contribution/DotPlot outputs) on a tiny per-cell-type subsample
# (<=40 cells/group) so API/argument errors surface in minutes instead
# of after a multi-hour full-scale run. Not meant to produce
# biologically meaningful results.

suppressMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(circlize)
  library(future)
  library(CellChat)
})
options(future.globals.maxSize = 1e12)
options(future.rng.onMisuse = "ignore")
set.seed(1)
plan("multicore", workers = 4)

Indir  <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/RenalUrotheliumScripts/output/"
Outdir <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/CellCommunicationScripts/output/smoketest/"
dir.create(Outdir, recursive = TRUE, showWarnings = FALSE)

CellTypeCol <- "final_annotation_with_uro"
UroLabel <- "Urothelium"

cat("Loading object...\n")
obj <- readRDS(paste0(Indir, "RenalUrothelium_allcells_scvi_annotations_meta_uro_updated.rds"))
DefaultAssay(obj) <- "RNA"

obj <- subset(obj, cells = colnames(obj)[!is.na(obj@meta.data[[CellTypeCol]])])
Idents(obj) <- obj@meta.data[[CellTypeCol]]
Idents(obj) <- droplevels(Idents(obj))

cat("Subsampling to <=40 cells per group for smoke test...\n")
set.seed(1)
cells.use <- unlist(lapply(split(colnames(obj), Idents(obj)), function(x) sample(x, min(40, length(x)))))
sub <- subset(obj, cells = cells.use)
Idents(sub) <- droplevels(Idents(sub))
print(table(Idents(sub)))

cat("createCellChat...\n")
cellchat <- createCellChat(object = sub, group.by = CellTypeCol, assay = "RNA")
cellchat@idents <- droplevels(cellchat@idents)

CellChatDB <- CellChatDB.mouse
cellchat@DB <- CellChatDB

cat("subsetData...\n")
cellchat <- subsetData(cellchat)
cat("identifyOverExpressedGenes...\n")
cellchat <- identifyOverExpressedGenes(cellchat)
cat("identifyOverExpressedInteractions...\n")
cellchat <- identifyOverExpressedInteractions(cellchat)
cat("computeCommunProb...\n")
cellchat <- computeCommunProb(cellchat, type = "triMean")
cat("filterCommunication...\n")
cellchat <- filterCommunication(cellchat, min.cells = 1)
cat("computeCommunProbPathway...\n")
cellchat <- computeCommunProbPathway(cellchat)
cat("aggregateNet...\n")
cellchat <- aggregateNet(cellchat)
cat("netAnalysis_computeCentrality...\n")
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

df.net <- subsetCommunication(cellchat)
cat("nrow df.net:", nrow(df.net), "\n")

groupSize <- as.numeric(table(cellchat@idents))
pdf(paste0(Outdir, "smoketest_circle.pdf"), height = 6, width = 12)
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(cellchat@net$count,  vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "count")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "weight")
dev.off()

pdf(paste0(Outdir, "smoketest_heatmap.pdf"), height = 8, width = 9)
print(netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues"))
dev.off()

AllPathways <- cellchat@netP$pathways
cat("Number of pathways detected:", length(AllPathways), "\n")

if (length(AllPathways) > 0) {
  ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", signaling = AllPathways, height = 10)
  ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", signaling = AllPathways, height = 10)
  pdf(paste0(Outdir, "smoketest_role_heatmap.pdf"), height = 10, width = 14)
  ht1 + ht2
  dev.off()
}

AllCellTypes <- levels(cellchat@idents)
OtherCellTypes <- setdiff(AllCellTypes, UroLabel)
UroAsSender <- df.net %>% filter(source == UroLabel)
UroAsReceiver <- df.net %>% filter(target == UroLabel)
cat("UroAsSender rows:", nrow(UroAsSender), " UroAsReceiver rows:", nrow(UroAsReceiver), "\n")

if (nrow(UroAsSender) > 0) {
  pdf(paste0(Outdir, "smoketest_uro_sender_bubble.pdf"), height = 6, width = 12)
  netVisual_bubble(cellchat, sources.use = UroLabel, targets.use = OtherCellTypes, remove.isolate = TRUE)
  dev.off()

  UroPathwaysAsSender <- unique(UroAsSender$pathway_name)
  pdf(paste0(Outdir, "smoketest_uro_sender_contrib.pdf"), height = 5, width = 6)
  ContribSender <- netAnalysis_contribution(cellchat, signaling = UroPathwaysAsSender, sources.use = UroLabel, targets.use = OtherCellTypes, return.data = TRUE)
  dev.off()
  cat("ContribSender LR.contribution rows:", nrow(ContribSender$LR.contribution), "\n")

  SenderLigands <- unique(UroAsSender$ligand)
  UroObj <- subset(sub, idents = UroLabel)
  p <- DotPlot(UroObj, features = rev(SenderLigands), cols = c("lightgrey", "blue")) + coord_flip()
  cat("DotPlot for sender ligands built OK, class:", class(p)[1], "\n")
}

cat("SMOKETEST COMPLETE - ALL STEPS RAN\n")
