################################################################################
# 10_UrotheliumDevelopment_chromVAR.R
#
# chromVAR TF motif accessibility, chromatin-level validation of the
# 06_UrotheliumDevelopment_SCENIC.R RNA-regulon story, for a fixed 20-TF list
# supplied by the user (Tead1, Esrrg, Atf2, E2f3, Tead2, Ezh2, Egr1, Rcor1,
# Jun, Fos, Runx1, Irf2, Klf7, Fli1, Xrcc4, Maf, Nfkb1, Nr1h4, Sp4, Stat1).
#
# JASPAR2020 species scope matters a lot here (checked interactively before
# writing this script): querying opts=list(species="Mus musculus") only
# recovers 4/20 of these TFs (107 mouse-only PFMs total). JASPAR motifs are
# cross-species-applicable sequence models, not organism-specific -- the
# standard-practice query is opts=list(tax_group="vertebrates",
# collection="CORE") (746 PFMs), which recovers 16/20. The remaining 4 --
# Ezh2 (Polycomb histone methyltransferase), Rcor1 (co-repressor/scaffold),
# Klf7, Xrcc4 (DNA repair factor) -- have NO JASPAR CORE vertebrate entry at
# all, confirmed by direct query, because none of them are sequence-specific
# DNA-binding factors (Klf7 is the one surprise -- the other 3 were expected
# on biological grounds). These 4 are reported as excluded, not silently
# dropped or forced.
#
# chromVAR's background-peak model (getBackgroundPeaks, inside RunChromVAR)
# samples ~50 GC-/accessibility-matched peaks from the FULL candidate peak
# set to build each peak's null distribution -- this is why 08 called peaks
# genome-wide-in-coordinate-space rather than restricting to loci near genes
# of interest; a peak set narrowed to a handful of loci would badly bias
# that background sampling.
#
# Row order for the heatmap: "early-shutoff on top, late-activating on
# bottom" is a monotonic ordering (by stage of peak/trough activity), not a
# clustering criterion -- pheatmap's dendrogram (used in 07 for the GO-term
# gene heatmap) won't reliably produce it, so cluster_rows=FALSE with rows
# pre-sorted by which.max/which.min stage index instead (see below).
#
# Input:  UrotheliumDevelopmentScripts/output/UrotheliumOnly_RNA_ATAC_WNN.rds
#         (written by 09_UrotheliumDevelopment_WNN.R)
# Output: UrotheliumDevelopment_chromVAR_TFcoverage.csv
#         UrotheliumDevelopment_chromVAR_zscores_by_stage.csv
#         FigC_UrotheliumDevelopment_chromVAR_TF_Heatmap.pdf
#         UrotheliumOnly_RNA_ATAC_chromVAR.rds
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(chromVAR)
  library(motifmatchr)
  library(TFBSTools)
  library(JASPAR2020)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
})

SCRIPT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/UrotheliumDevelopmentScripts"
OUT_DIR    <- file.path(SCRIPT_DIR, "output")

STAGE_ORDER <- c("E16.5", "P0", "W3", "W12", "W52", "W92")
TF_PANEL <- c("Tead1", "Esrrg", "Atf2", "E2f3", "Tead2", "Ezh2", "Egr1", "Rcor1",
              "Jun", "Fos", "Runx1", "Irf2", "Klf7", "Fli1", "Xrcc4", "Maf",
              "Nfkb1", "Nr1h4", "Sp4", "Stat1")

# ── Load ─────────────────────────────────────────────────────────────────────
message("==> Loading combined RNA+ATAC+WNN Uro object ...")
in_rds <- file.path(OUT_DIR, "UrotheliumOnly_RNA_ATAC_WNN.rds")
if (!file.exists(in_rds)) {
  stop("Missing ", in_rds, " -- run 09_UrotheliumDevelopment_WNN.R first.")
}
uro <- readRDS(in_rds)
uro$Age <- factor(uro$Age, levels = STAGE_ORDER)
DefaultAssay(uro) <- "ATAC"

################################################################################
# JASPAR motifs: vertebrate CORE, not mouse-only (see header)
################################################################################
message("\n==> Fetching JASPAR2020 vertebrate CORE PFMs ...")
pfm <- getMatrixSet(JASPAR2020, list(tax_group = "vertebrates", collection = "CORE"))
message(sprintf("  %d PFMs fetched", length(pfm)))

pfm_names <- sapply(pfm, function(x) name(x))
tf_hit <- sapply(TF_PANEL, function(g) any(grepl(paste0("(^|::)", g, "(::|$)"), pfm_names, ignore.case = TRUE)))
coverage_df <- data.frame(TF = TF_PANEL, has_JASPAR_motif = tf_hit)
write.csv(coverage_df, file.path(OUT_DIR, "UrotheliumDevelopment_chromVAR_TFcoverage.csv"), row.names = FALSE)
message(sprintf("  %d / %d requested TFs have a JASPAR vertebrate CORE motif", sum(tf_hit), length(TF_PANEL)))
message("  Excluded (no JASPAR motif -- not sequence-specific DNA-binding factors): ",
        paste(TF_PANEL[!tf_hit], collapse = ", "))
TF_PANEL_MATCHED <- TF_PANEL[tf_hit]

################################################################################
# AddMotifs + RunChromVAR
################################################################################
message("\n==> Scanning peaks for motif matches (AddMotifs) ...")
uro <- AddMotifs(uro, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pfm, assay = "ATAC")

message("\n==> Running chromVAR (RunChromVAR) ...")
uro <- RunChromVAR(uro, genome = BSgenome.Mmusculus.UCSC.mm10, assay = "ATAC")
DefaultAssay(uro) <- "chromvar"

message(sprintf("  chromVAR assay: %d motifs x %d cells", nrow(uro[["chromvar"]]), ncol(uro[["chromvar"]])))

################################################################################
# Per-stage mean z-score for the matched TF panel
################################################################################
message("\n==> Computing per-stage mean chromVAR z-score for matched TFs ...")

motif_meta <- uro[["ATAC"]]@motifs@motif.names
motif_ids_for_tf <- sapply(TF_PANEL_MATCHED, function(g) {
  hit <- names(motif_meta)[sapply(motif_meta, function(nm) grepl(paste0("(^|::)", g, "(::|$)"), nm, ignore.case = TRUE))]
  if (length(hit) == 0) return(NA_character_)
  hit[1]  # first match if a TF has multiple JASPAR entries
})
motif_ids_for_tf <- motif_ids_for_tf[!is.na(motif_ids_for_tf)]
message(sprintf("  %d / %d matched TFs resolved to a motif ID present in this peak set's scan",
                 length(motif_ids_for_tf), length(TF_PANEL_MATCHED)))

zmat_cells <- GetAssayData(uro, assay = "chromvar", layer = "data")[motif_ids_for_tf, , drop = FALSE]
rownames(zmat_cells) <- names(motif_ids_for_tf)

zmat_stage <- sapply(STAGE_ORDER, function(st) {
  Matrix::rowMeans(zmat_cells[, uro$Age == st, drop = FALSE], na.rm = TRUE)
})
rownames(zmat_stage) <- rownames(zmat_cells)
colnames(zmat_stage) <- STAGE_ORDER

write.csv(as.data.frame(zmat_stage) %>% tibble::rownames_to_column("TF"),
          file.path(OUT_DIR, "UrotheliumDevelopment_chromVAR_zscores_by_stage.csv"), row.names = FALSE)
message("  Saved: UrotheliumDevelopment_chromVAR_zscores_by_stage.csv")

saveRDS(uro, file.path(OUT_DIR, "UrotheliumOnly_RNA_ATAC_chromVAR.rds"))
message("  Saved: UrotheliumOnly_RNA_ATAC_chromVAR.rds")

################################################################################
# Figure C: TF motif activity heatmap
################################################################################
message("\n==> Building Figure C (chromVAR TF heatmap) ...")

# Row-scale (z-score across stages, per TF) then order rows by which stage
# each TF peaks (early-shutoff on top, late-activating on bottom) -- a
# monotonic criterion, not a clustering one (see header).
zmat_rowscaled <- t(scale(t(zmat_stage)))
peak_stage_idx <- apply(zmat_rowscaled, 1, which.max)
row_order <- order(peak_stage_idx)
zmat_ordered <- zmat_stage[row_order, , drop = FALSE]

hm_palette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
hm_breaks  <- seq(-2.5, 2.5, length.out = length(hm_palette) + 1)

pdf(file.path(OUT_DIR, "FigC_UrotheliumDevelopment_chromVAR_TF_Heatmap.pdf"),
    width = 6, height = max(4, 0.3 * nrow(zmat_ordered)))
pheatmap(zmat_ordered, scale = "row", cluster_rows = FALSE, cluster_cols = FALSE,
         color = hm_palette, breaks = hm_breaks, border_color = NA,
         main = sprintf("Panel C. chromVAR TF motif activity by stage (%d/%d TFs matched)",
                         nrow(zmat_ordered), length(TF_PANEL)))
dev.off()
message("  Saved: FigC_UrotheliumDevelopment_chromVAR_TF_Heatmap.pdf")

message("\n==> Done.")
