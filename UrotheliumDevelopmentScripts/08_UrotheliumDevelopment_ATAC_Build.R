################################################################################
# 08_UrotheliumDevelopment_ATAC_Build.R
#
# Attach the matched multiome ATAC modality to the 715-cell Uro RNA object,
# call peaks, and build the combined RNA+ATAC object that scripts 09-11
# (WNN, chromVAR, LinkPeaks) all load from.
#
# The Chen2025 RDS/h5ad on this filesystem is RNA-only despite the
# "MultiOmic" name -- CELLxGENE only hosts the RNA half of multiome deposits.
# The user separately supplied the raw ATAC fragments file. Confirmed
# interactively before writing this script:
#   - Genome build is mm10 (paper's own Methods: CellRanger-ARC v2.0.0,
#     mouse reference mm10).
#   - Fragments file barcodes (e.g. "20210319-NMKE16_5F_AAACAGCCAAACAACA-1")
#     exactly match colnames(uro) (sample_id + "_" + 16bp barcode + "-1") --
#     no remapping needed.
#   - Chromosome naming is UCSC-style ("chr1"...), 21 seqnames, no chrM.
#   - No peaks/WNN/chromVAR/anything ATAC-derived exists anywhere already.
#
# CallPeaks() gotcha (caught in a Plan-agent review of an earlier draft,
# confirmed against Signac source): CallPeaks() does NOT respect a Fragment
# object's cells= whitelist -- CallPeaks.ChromatinAssay pulls the raw
# fragment file PATH and shells out to MACS on the whole file regardless of
# what's attached to the Seurat object. The only cell-aware path is
# CallPeaks.Seurat(object, group.by=...), which uses SplitFragments()
# internally to stream the file once and write per-group filtered fragment
# files before invoking MACS per group. Calling group.by = "Age" here (not
# a single flat "Uro" group) additionally avoids MACS's default human
# effective-genome-size (2.7e9) -- overridden to "mm" (1.87e9) below -- and
# avoids diluting this rare population's peaks by pooling with the other
# ~202,400 non-Uro atlas cells, which is what a genome-wide-on-the-full-atlas
# peak call would have done. FeatureMatrix()'s cells= argument DOES work as
# documented (restricts output columns), unlike CallPeaks()'s.
#
# EnsDb.Mmusculus.v79 uses Ensembl-style seqnames (1, 2, ..., X) while the
# fragments file/peaks use UCSC-style ("chr1", ...) -- seqlevelsStyle is
# forced to "UCSC" before attaching gene annotations, or Annotation()
# assignment silently produces an empty/mismatched gene track and 11's
# CoveragePlot gene model would just be blank.
#
# RegionStats (GC content / width per peak) is computed once here, not
# repeated in 10/chromVAR or 11/LinkPeaks, since both need it and it's
# deterministic given the peak set.
#
# Input:  UrotheliumDevelopmentScripts/output/UrotheliumOnly_full_annotated.rds
#         (written by 05_UrotheliumDevelopment_ModuleScores.R; carries the
#         Uro-specific "pca" reduction from 02 and seurat_clusters used by 09)
#         /home/gdbecknelllab/xxw004/gdjacksonlab/RenalDiseases/Mouse/
#           DevelopDefect/MultiOmic_KidneyDevNatGen/
#           MultiOmicSpatialMouseKidney_SexSpecific_lifespan_Chen2025NatGenet_ATAC.tsv.bgz
# Output: UrotheliumOnly_RNA_ATAC.rds
#         UrotheliumDevelopment_ATAC_PeaksPerStage.csv
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(EnsDb.Mmusculus.v79)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(dplyr)
})

SCRIPT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/UrotheliumDevelopmentScripts"
OUT_DIR    <- file.path(SCRIPT_DIR, "output")

STAGE_ORDER <- c("E16.5", "P0", "W3", "W12", "W52", "W92")

FRAGMENTS_PATH <- "/home/gdbecknelllab/xxw004/gdjacksonlab/RenalDiseases/Mouse/DevelopDefect/MultiOmic_KidneyDevNatGen/MultiOmicSpatialMouseKidney_SexSpecific_lifespan_Chen2025NatGenet_ATAC.tsv.bgz"

# Resolved via absolute path (MACS3_BIN env var set by the submit script),
# NOT PATH lookup -- scrna_env (where macs3 lives) ships its own R/Rscript
# that would silently shadow this module-loaded R/4.4.0 if its bin/ were
# ever prepended to PATH, running the whole job under the wrong R install
# (caught interactively before this ever reached a submitted job).
macs3_bin <- Sys.getenv("MACS3_BIN", unset = Sys.which("macs3"))
if (macs3_bin == "") {
  stop("macs3 not found -- submit script must set MACS3_BIN to its absolute path.")
}
message("Using MACS3 at: ", macs3_bin)

# ── Load the canonical Uro RNA object ───────────────────────────────────────
message("==> Loading Uro RNA object ...")
in_rds <- file.path(OUT_DIR, "UrotheliumOnly_full_annotated.rds")
if (!file.exists(in_rds)) {
  stop("Missing ", in_rds, " -- run 05_UrotheliumDevelopment_ModuleScores.R first.")
}
uro <- readRDS(in_rds)
uro$Age <- factor(uro$Age, levels = STAGE_ORDER)
DefaultAssay(uro) <- "originalexp"
message(sprintf("  %d Uro cells, reductions: %s", ncol(uro), paste(Reductions(uro), collapse = ", ")))

# ── Attach ATAC fragments (no peaks/counts yet) ─────────────────────────────
message("\n==> Attaching ATAC fragments ...")
if (!file.exists(FRAGMENTS_PATH)) {
  stop("Missing ATAC fragments file: ", FRAGMENTS_PATH)
}
uro_frag <- CreateFragmentObject(path = FRAGMENTS_PATH, cells = colnames(uro))

# CreateChromatinAssay requires *some* counts matrix with valid
# "chr-start-end" rownames to build its internal ranges -- this placeholder
# single-region matrix exists only so CallPeaks(group.by="Age") below has an
# assay to find the attached Fragment object via Fragments(object[[assay]]);
# it gets replaced wholesale by the real peak x cell matrix right after.
placeholder_region <- "chr1-3000000-3000001"
atac_placeholder <- CreateChromatinAssay(
  counts   = Matrix::Matrix(0, nrow = 1, ncol = ncol(uro),
                             dimnames = list(placeholder_region, colnames(uro)), sparse = TRUE),
  fragments = list(uro_frag),
  min.cells = 0, min.features = 0
)
uro[["ATAC"]] <- atac_placeholder
DefaultAssay(uro) <- "ATAC"

# ── Call peaks: cell-aware path only, per-stage then merged ────────────────
message("\n==> Calling peaks (MACS3, group.by = Age, mouse effective genome size) ...")
peaks <- CallPeaks(
  object                 = uro,
  assay                  = "ATAC",
  group.by               = "Age",
  macs2.path             = macs3_bin,
  effective.genome.size  = 1.87e9,   # "mm" -- MACS's default (2.7e9) is human
  combine.peaks           = TRUE,
  outdir                 = file.path(OUT_DIR, "macs3_tmp"),
  cleanup                = TRUE
)
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subset(peaks, seqnames %in% paste0("chr", c(1:19, "X", "Y")))

message(sprintf("  %d merged peaks across %d stage(s)", length(peaks), length(STAGE_ORDER)))
if ("peak_called_in" %in% colnames(mcols(peaks))) {
  per_stage_counts <- sapply(STAGE_ORDER, function(st) sum(grepl(st, peaks$peak_called_in, fixed = TRUE)))
  message("  Peaks called in each stage (non-exclusive, a peak can be called in >1 stage):")
  print(per_stage_counts)
  write.csv(data.frame(stage = STAGE_ORDER, n_peaks_called = per_stage_counts),
            file.path(OUT_DIR, "UrotheliumDevelopment_ATAC_PeaksPerStage.csv"), row.names = FALSE)
}

# ── Real counts: FeatureMatrix, Uro cells only ──────────────────────────────
message("\n==> Building peak x cell count matrix (Uro cells only) ...")
counts <- FeatureMatrix(fragments = Fragments(uro[["ATAC"]]), features = peaks, cells = colnames(uro))

atac_assay <- CreateChromatinAssay(
  counts    = counts,
  fragments = list(uro_frag),
  min.cells = 0, min.features = 0
)
uro[["ATAC"]] <- atac_assay
DefaultAssay(uro) <- "ATAC"
message(sprintf("  ATAC assay: %d peaks x %d cells", nrow(uro[["ATAC"]]), ncol(uro[["ATAC"]])))

# ── Gene annotation (for CoveragePlot gene models in script 11) ────────────
message("\n==> Attaching gene annotation (EnsDb.Mmusculus.v79, UCSC seqlevels) ...")
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "mm10"
Annotation(uro[["ATAC"]]) <- annotations

# ── RegionStats: GC content / width per peak, needed by 10 and 11 ─────────
message("\n==> Computing RegionStats (GC content, width) ...")
uro <- RegionStats(uro, genome = BSgenome.Mmusculus.UCSC.mm10, assay = "ATAC")

# ── QC ───────────────────────────────────────────────────────────────────
message("\n==> QC ...")
uro <- TSSEnrichment(uro, assay = "ATAC", fast = TRUE)
message("  TSS enrichment summary (typical good ATAC data: ~5-15):")
print(summary(uro$TSS.enrichment))
message("  Peak counts per cell x stage summary:")
uro$atac_total_counts <- colSums(GetAssayData(uro, assay = "ATAC", layer = "counts"))
print(uro@meta.data %>% group_by(Age) %>%
        summarise(median_atac_counts = median(atac_total_counts), n = n(), .groups = "drop"))

# ── Save ─────────────────────────────────────────────────────────────────
saveRDS(uro, file.path(OUT_DIR, "UrotheliumOnly_RNA_ATAC.rds"))
message("\n  Saved: UrotheliumOnly_RNA_ATAC.rds")

message("\n==> Done.")
