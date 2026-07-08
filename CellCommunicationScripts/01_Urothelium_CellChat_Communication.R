#!/usr/bin/env Rscript
# Cell-cell communication: Urothelium to other kidney cell types
# Author: Xin Wang
# Email: xin.wang@nationwidechildrens.org
# Copyright (c) 2026 Kidney and Urology Tract Center, Nationwide Children's Hospital
#
# Description:
#   Profiles cell-cell communication across all annotated kidney cell types in
#   RenalUrothelium_allcells_scvi_annotations_meta_uro_updated.rds
#   (cell type column: final_annotation_with_uro), restricted to healthy-kidney
#   cells only (FinalConditionL1 == "HealthyKidney"), then extracts and
#   visualizes the Urothelium <-> other-cell-type interactome specifically.
#
# Method: CellChat (https://github.com/jinworks/CellChat)
# Pipeline structure follows:
#   SingleCell_Day0_CellCommunication_V02062026.Rmd
#   Developed_GroupCompaision_CellCommunication_RequestedCellType_RequestedPathways_V02052026.Rmd
#
# Note: run with no downsampling beyond the HealthyKidney condition filter (per user request).

suppressMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(circlize)
  library(future)
  library(CellChat)
})

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 1e12)
options(future.rng.onMisuse = "ignore")
set.seed(200000)

# multicore (fork-based) avoids serializing the ~1M-cell object across workers,
# which multisession would otherwise require
n_workers <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))
plan("multicore", workers = n_workers)

# ############################################################
# 0. Paths and parameters
# ############################################################
Indir  <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/RenalUrotheliumScripts/output/"
Outdir <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/CellCommunicationScripts/output/"
UroDir <- paste0(Outdir, "Urothelium/")
dir.create(Outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(UroDir, recursive = TRUE, showWarnings = FALSE)
setwd(Outdir)

CellTypeCol   <- "final_annotation_with_uro"
UroLabel      <- "Urothelium"
ConditionCol  <- "FinalConditionL1"
KeepCondition <- "HealthyKidney"

# ############################################################
# 1. Read the annotated single-cell object
# ############################################################
SingleCellObj <- readRDS(paste0(Indir, "RenalUrothelium_allcells_scvi_annotations_meta_uro_updated.rds"))
DefaultAssay(SingleCellObj) <- "RNA"
if (length(Layers(SingleCellObj[["RNA"]])) > 2) {
  SingleCellObj <- JoinLayers(SingleCellObj)
}

# restrict to healthy-kidney cells only
SingleCellObj <- subset(SingleCellObj, cells = colnames(SingleCellObj)[SingleCellObj@meta.data[[ConditionCol]] == KeepCondition])
print(table(SingleCellObj@meta.data[[ConditionCol]]))

# drop any cells missing the cell type annotation
KeepCells <- colnames(SingleCellObj)[!is.na(SingleCellObj@meta.data[[CellTypeCol]])]
SingleCellObj <- subset(SingleCellObj, cells = KeepCells)
Idents(SingleCellObj) <- SingleCellObj@meta.data[[CellTypeCol]]
Idents(SingleCellObj) <- droplevels(Idents(SingleCellObj))
print(table(Idents(SingleCellObj)))

# ############################################################
# 2. Create the CellChat object (all cell types)
# ############################################################
cellchat <- createCellChat(object = SingleCellObj, group.by = CellTypeCol, assay = "RNA")
cellchat@idents <- droplevels(cellchat@idents)
print(table(cellchat@idents))

CellChatDB <- CellChatDB.mouse
pdf(paste0(Outdir, "CellChatDB_Category_Overview.pdf"), height = 6, width = 8)
showDatabaseCategory(CellChatDB)
dev.off()

cellchat@DB <- CellChatDB

# ############################################################
# 3. Preprocessing
# ############################################################
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# ############################################################
# 4. Inference of cell-cell communication network
# ############################################################
cellchat <- computeCommunProb(cellchat, type = "triMean")
saveRDS(cellchat, file = paste0(Outdir, "CellChat_Urothelium_AllCellTypes_raw.rds"))

cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

saveRDS(cellchat, file = paste0(Outdir, "CellChat_Urothelium_AllCellTypes_final.rds"))

df.net <- subsetCommunication(cellchat)
write.csv(df.net, file = paste0(Outdir, "CellChat_AllCellTypes_LR_interactions.csv"), row.names = FALSE)

# ############################################################
# 5. Global overview (all cell types)
# ############################################################
groupSize <- as.numeric(table(cellchat@idents))

pdf(paste0(Outdir, "Global_Interactions_CountWeight_Circle.pdf"), height = 6, width = 12)
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(cellchat@net$count,  vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights/strength")
dev.off()

pdf(paste0(Outdir, "Global_Interactions_Heatmap.pdf"), height = 9, width = 10)
print(netVisual_heatmap(cellchat, measure = "count",  color.heatmap = "Blues"))
print(netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Blues"))
dev.off()

AllPathways <- cellchat@netP$pathways
write.csv(data.frame(pathway = AllPathways), file = paste0(Outdir, "Global_DetectedPathways.csv"), row.names = FALSE)

ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", signaling = AllPathways, height = 20)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", signaling = AllPathways, height = 20)
pdf(paste0(Outdir, "Global_SignalingRole_Outgoing_Incoming_Heatmap.pdf"), height = 20, width = 18)
ht1 + ht2
dev.off()

pdf(paste0(Outdir, "Global_SignalingRole_Scatter.pdf"), height = 6, width = 7)
netAnalysis_signalingRole_scatter(cellchat)
dev.off()

# ############################################################
# 6. Urothelium-focused communication
# ############################################################
AllCellTypes   <- levels(cellchat@idents)
OtherCellTypes <- setdiff(AllCellTypes, UroLabel)

# 6a. Full LR interaction tables: Urothelium as sender / receiver
UroAsSender   <- df.net %>% filter(source == UroLabel)
UroAsReceiver <- df.net %>% filter(target == UroLabel)
write.csv(UroAsSender,   file = paste0(UroDir, "Urothelium_AsSender_LR_interactions.csv"),   row.names = FALSE)
write.csv(UroAsReceiver, file = paste0(UroDir, "Urothelium_AsReceiver_LR_interactions.csv"), row.names = FALSE)

# 6b. Circle plots isolating the Urothelium row / column
matWeight <- cellchat@net$weight

matSender <- matrix(0, nrow = nrow(matWeight), ncol = ncol(matWeight), dimnames = dimnames(matWeight))
matSender[UroLabel, ] <- matWeight[UroLabel, ]
matReceiver <- matrix(0, nrow = nrow(matWeight), ncol = ncol(matWeight), dimnames = dimnames(matWeight))
matReceiver[, UroLabel] <- matWeight[, UroLabel]

pdf(paste0(UroDir, "Urothelium_Sender_Receiver_Circle.pdf"), height = 6, width = 12)
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(matSender,   vertex.weight = groupSize, weight.scale = TRUE, edge.weight.max = max(matWeight), title.name = "Urothelium as sender")
netVisual_circle(matReceiver, vertex.weight = groupSize, weight.scale = TRUE, edge.weight.max = max(matWeight), title.name = "Urothelium as receiver")
dev.off()

# 6c. Bubble plots: Urothelium -> all other cell types, and all other cell types -> Urothelium
nSenderPairs   <- length(unique(UroAsSender$interaction_name_2))
nReceiverPairs <- length(unique(UroAsReceiver$interaction_name_2))

pdf(paste0(UroDir, "Urothelium_AsSender_Bubble.pdf"), height = max(6, nSenderPairs / 4), width = 12)
netVisual_bubble(cellchat, sources.use = UroLabel, targets.use = OtherCellTypes, remove.isolate = TRUE)
dev.off()

pdf(paste0(UroDir, "Urothelium_AsReceiver_Bubble.pdf"), height = max(6, nReceiverPairs / 4), width = 12)
netVisual_bubble(cellchat, sources.use = OtherCellTypes, targets.use = UroLabel, remove.isolate = TRUE)
dev.off()

# 6d. Pathway / ligand-receptor contribution driving Urothelium signaling
UroPathwaysAsSender   <- unique(UroAsSender$pathway_name)
UroPathwaysAsReceiver <- unique(UroAsReceiver$pathway_name)

pdf(paste0(UroDir, "Urothelium_AsSender_PathwayContribution.pdf"), height = max(5, length(UroPathwaysAsSender) / 3), width = 6)
ContribSender <- netAnalysis_contribution(cellchat, signaling = UroPathwaysAsSender, sources.use = UroLabel, targets.use = OtherCellTypes, return.data = TRUE)
dev.off()
write.csv(ContribSender$LR.contribution, file = paste0(UroDir, "Urothelium_AsSender_LR_Contribution.csv"), row.names = FALSE)

pdf(paste0(UroDir, "Urothelium_AsReceiver_PathwayContribution.pdf"), height = max(5, length(UroPathwaysAsReceiver) / 3), width = 6)
ContribReceiver <- netAnalysis_contribution(cellchat, signaling = UroPathwaysAsReceiver, sources.use = OtherCellTypes, targets.use = UroLabel, return.data = TRUE)
dev.off()
write.csv(ContribReceiver$LR.contribution, file = paste0(UroDir, "Urothelium_AsReceiver_LR_Contribution.csv"), row.names = FALSE)

# 6e. Total communication strength between Urothelium and each other cell type
UroStrengthSender <- UroAsSender %>%
  group_by(target) %>%
  summarise(total_prob = sum(prob)) %>%
  arrange(desc(total_prob))
write.csv(UroStrengthSender, file = paste0(UroDir, "Urothelium_AsSender_StrengthByTarget.csv"), row.names = FALSE)

UroBarSender <- ggplot(UroStrengthSender, aes(x = reorder(target, total_prob), y = total_prob)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() + theme_classic() +
  labs(x = "Target cell type", y = "Total communication probability", title = "Urothelium -> other cell types")
ggsave(paste0(UroDir, "Urothelium_AsSender_StrengthByTarget_Barplot.pdf"), plot = UroBarSender, height = 8, width = 6)

UroStrengthReceiver <- UroAsReceiver %>%
  group_by(source) %>%
  summarise(total_prob = sum(prob)) %>%
  arrange(desc(total_prob))
write.csv(UroStrengthReceiver, file = paste0(UroDir, "Urothelium_AsReceiver_StrengthBySource.csv"), row.names = FALSE)

UroBarReceiver <- ggplot(UroStrengthReceiver, aes(x = reorder(source, total_prob), y = total_prob)) +
  geom_bar(stat = "identity", fill = "firebrick") +
  coord_flip() + theme_classic() +
  labs(x = "Source cell type", y = "Total communication probability", title = "Other cell types -> Urothelium")
ggsave(paste0(UroDir, "Urothelium_AsReceiver_StrengthBySource_Barplot.pdf"), plot = UroBarReceiver, height = 8, width = 6)

# 6f. Ligand / receptor gene expression DotPlots
UroObj   <- subset(SingleCellObj, idents = UroLabel)
OtherObj <- subset(SingleCellObj, idents = OtherCellTypes)
Idents(OtherObj) <- factor(Idents(OtherObj), levels = OtherCellTypes)

## Urothelium as sender: ligand genes in Urothelium, receptor genes in target cell types
SenderLigands   <- unique(UroAsSender$ligand)
SenderReceptors <- unique(UroAsSender$receptor)

UroLigandDotPlot <- DotPlot(UroObj, features = rev(SenderLigands), cols = c("lightgrey", "blue")) +
  coord_flip() + theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8), axis.title = element_blank())
OtherReceptorDotPlot <- DotPlot(OtherObj, features = rev(SenderReceptors), cols = c("lightgrey", "red")) +
  coord_flip() + theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8), axis.title = element_blank())

ggsave(paste0(UroDir, "Urothelium_AsSender_Ligand_DotPlot.pdf"), plot = UroLigandDotPlot, height = max(4, length(SenderLigands) / 4), width = 5)
ggsave(paste0(UroDir, "Urothelium_AsSender_TargetReceptor_DotPlot.pdf"), plot = OtherReceptorDotPlot, height = max(4, length(SenderReceptors) / 4), width = 10)

## Urothelium as receiver: ligand genes in other cell types, receptor genes in Urothelium
ReceiverLigands   <- unique(UroAsReceiver$ligand)
ReceiverReceptors <- unique(UroAsReceiver$receptor)

OtherLigandDotPlot <- DotPlot(OtherObj, features = rev(ReceiverLigands), cols = c("lightgrey", "blue")) +
  coord_flip() + theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8), axis.title = element_blank())
UroReceptorDotPlot <- DotPlot(UroObj, features = rev(ReceiverReceptors), cols = c("lightgrey", "red")) +
  coord_flip() + theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8), axis.title = element_blank())

ggsave(paste0(UroDir, "Urothelium_AsReceiver_SourceLigand_DotPlot.pdf"), plot = OtherLigandDotPlot, height = max(4, length(ReceiverLigands) / 4), width = 10)
ggsave(paste0(UroDir, "Urothelium_AsReceiver_Receptor_DotPlot.pdf"), plot = UroReceptorDotPlot, height = max(4, length(ReceiverReceptors) / 4), width = 5)

sessionInfo()
