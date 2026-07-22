################################################################################
# 13_UrotheliumDevelopment_GenePanel_BigWig.R
#
# Per-stage, genome-wide BigWig tracks (Uro cells only), for loading into an
# external genome browser (IGV/UCSC) alongside the static CoveragePlot pages
# from 12_UrotheliumDevelopment_GenePanel_CoveragePlots.R.
#
# Signac's SplitFragments() streams the fragments file once, writing one
# filtered-and-grouped file per Age level -- the fragments file it reads
# from is already restricted to these 715 Uro cells (CreateFragmentObject's
# cells= whitelist in 08), so "per stage" here already means "per stage,
# Uro cells only" without any further cell filtering needed.
#
# Normalization: each stage's track is scaled to fragments-per-million
# (total valid fragments in that stage's split file), the standard
# convention for cross-sample-comparable ATAC bigwigs -- not the same
# normalization CoverageTrack uses internally for 12's plots (that one is
# also group-size- and smoothing-window-dependent), so track heights here
# and in 12's PDF are each internally consistent but not numerically
# identical to one another.
#
# Input:  UrotheliumDevelopmentScripts/output/UrotheliumOnly_RNA_ATAC_chromVAR.rds
#         (written by 10_UrotheliumDevelopment_chromVAR.R; carries the
#         Fragment object originally attached in 08, whitelisted to these
#         715 Uro cells)
# Output: UrotheliumDevelopment_ATAC_<stage>.bw   (one per stage, 6 total,
#         genome-wide coverage)
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(rtracklayer)
  library(data.table)
  library(dplyr)
})

SCRIPT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/UrotheliumDevelopmentScripts"
OUT_DIR    <- file.path(SCRIPT_DIR, "output")
SPLIT_DIR  <- file.path(OUT_DIR, "bigwig_split_fragments")
dir.create(SPLIT_DIR, showWarnings = FALSE, recursive = TRUE)

STAGE_ORDER <- c("E16.5", "P0", "W3", "W12", "W52", "W92")

# ── Load ─────────────────────────────────────────────────────────────────────
message("==> Loading combined RNA+ATAC+chromVAR Uro object ...")
in_rds <- file.path(OUT_DIR, "UrotheliumOnly_RNA_ATAC_chromVAR.rds")
if (!file.exists(in_rds)) {
  stop("Missing ", in_rds, " -- run 10_UrotheliumDevelopment_chromVAR.R first.")
}
uro <- readRDS(in_rds)
uro$Age <- factor(uro$Age, levels = STAGE_ORDER)
DefaultAssay(uro) <- "ATAC"

# ── Split fragments by stage (Uro cells only, whole file scanned once) ─────
message("\n==> Splitting fragments by stage (Uro cells only) ...")
SplitFragments(
  object   = uro,
  assay    = "ATAC",
  group.by = "Age",
  outdir   = SPLIT_DIR,
  verbose  = TRUE
)
split_files <- list.files(SPLIT_DIR, pattern = "\\.bed$", full.names = TRUE)
message("  Split files written: ", paste(basename(split_files), collapse = ", "))

# ── Per-stage BigWig: genome-wide coverage ──────────────────────────────────
message("\n==> Building per-stage BigWig tracks (genome-wide) ...")
seqlengths_mm10 <- seqlengths(BSgenome.Mmusculus.UCSC.mm10)

for (st in STAGE_ORDER) {
  frag_file <- file.path(SPLIT_DIR, paste0(st, ".bed"))
  if (!file.exists(frag_file)) {
    message(sprintf("  [skip] %s: no split fragment file found", st))
    next
  }
  message(sprintf("  %s ...", st))
  frags <- fread(frag_file, header = FALSE)
  col_names <- c("chrom", "start", "end", "barcode", "count")
  setnames(frags, seq_len(min(ncol(frags), length(col_names))), col_names[seq_len(min(ncol(frags), length(col_names)))])

  gr <- GRanges(frags$chrom, IRanges(frags$start + 1, frags$end))
  total_fragments <- length(gr)

  chroms_present <- intersect(seqlevels(gr), names(seqlengths_mm10))
  gr <- keepSeqlevels(gr, chroms_present, pruning.mode = "coarse")
  seqlengths(gr) <- seqlengths_mm10[chroms_present]
  message(sprintf("    %d fragments, %d chromosomes", total_fragments, length(chroms_present)))

  cov <- coverage(gr)
  # Fragments-per-million normalization (standard cross-sample-comparable
  # convention for ATAC bigwigs).
  scale_factor <- 1e6 / total_fragments
  cov <- cov * scale_factor

  out_bw <- file.path(OUT_DIR, sprintf("UrotheliumDevelopment_ATAC_%s.bw", st))
  export.bw(cov, out_bw)
  message("    Saved: ", out_bw)
}

message("\n==> Done.")
