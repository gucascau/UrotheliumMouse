################################################################################
# 12_UrotheliumDevelopment_GenePanel_CoveragePlots.R
#
# ATAC coverage/peak-track plots (Signac::CoveragePlot), one page per gene,
# across the 6 developmental stages, for a 13-gene panel spanning Spp1
# (osteopontin, this pipeline's main CellChat/LinkPeaks focus), the
# "Upper-tract identity" module (Pax8, Pax2, Glis3, Fgfr2, Pkhd1, Bicc1 --
# same set as 01/05/06's module score), Col4a3/Col4a4 (basement-membrane
# collagen IV chains, see CellCommunicationScripts/02's Uro-only DotPlot),
# Trp53, and uroplakins/keratins (Upk1b, Upk3a, Krt20).
#
# Purely descriptive coverage tracks -- no LinkPeaks correlation/links layer
# here (unlike 11's Panel D), since LinkPeaks was only ever run against
# 11's own smaller gene panel and that run's result isn't persisted to disk
# (11 never calls saveRDS()); re-running LinkPeaks for this different
# 13-gene panel just to overlay links wasn't asked for.
#
# All 13 pages share one y-axis scale (CoveragePlot's ymax) rather than each
# auto-scaling to its own local max, so peak heights are visually comparable
# across genes -- see the "Shared y-axis scale" section below for how that
# value is computed (Signac's own CoverageTrack() normalization formula,
# reproduced from source, not a fixed/eyeballed number).
#
# expression.assay: originalexp's rownames are Ensembl IDs throughout this
# whole project (same recurring gotcha as 11_UrotheliumDevelopment_
# LinkPeaks.R and CellCommunicationScripts/02's Spp1/Col4a3-5 DotPlot) --
# CoveragePlot's features= argument needs symbol rownames to match this
# panel's gene symbols, and rownames() can't be reassigned on the object
# directly ("Renaming features in v3/v4 assays is not supported"), so a
# small symbol-keyed Assay is built for just this panel, same fix as 11.
#
# Input:  UrotheliumDevelopmentScripts/output/UrotheliumOnly_RNA_ATAC_chromVAR.rds
#         (written by 10_UrotheliumDevelopment_chromVAR.R)
# Output: FigG_UrotheliumDevelopment_GenePanel_CoveragePlots.pdf
#         (one page per gene, genes not found in the panel reported and
#         skipped rather than erroring the whole run)
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(ggplot2)
  library(patchwork)
  library(viridisLite)
  library(dplyr)
})

SCRIPT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/UrotheliumDevelopmentScripts"
OUT_DIR    <- file.path(SCRIPT_DIR, "output")

STAGE_ORDER  <- c("E16.5", "P0", "W3", "W12", "W52", "W92")
STAGE_COLORS <- setNames(viridisLite::viridis(length(STAGE_ORDER)), STAGE_ORDER)

GENE_PANEL <- unique(c(
  "Spp1",
  "Pax8", "Pax2", "Glis3", "Fgfr2", "Pkhd1", "Bicc1",
  "Col4a3", "Col4a4",
  "Trp53",
  "Upk1b", "Upk3a", "Krt20"
))

# ── Load ─────────────────────────────────────────────────────────────────────
message("==> Loading combined RNA+ATAC+chromVAR Uro object ...")
in_rds <- file.path(OUT_DIR, "UrotheliumOnly_RNA_ATAC_chromVAR.rds")
if (!file.exists(in_rds)) {
  stop("Missing ", in_rds, " -- run 10_UrotheliumDevelopment_chromVAR.R first.")
}
uro <- readRDS(in_rds)
uro$Age <- factor(uro$Age, levels = STAGE_ORDER)
DefaultAssay(uro) <- "ATAC"

fm_rna <- uro[["originalexp"]]@meta.features
annotation_genes <- unique(Annotation(uro[["ATAC"]])$gene_name)

in_expression <- GENE_PANEL %in% fm_rna$gene_symbols
in_annotation <- GENE_PANEL %in% annotation_genes
plottable     <- in_expression & in_annotation

message(sprintf("Genes with both RNA expression and ATAC gene-model coordinates: %d / %d",
                sum(plottable), length(GENE_PANEL)))
if (any(!in_expression)) {
  message("  Not found in RNA gene_symbols (skipped): ", paste(GENE_PANEL[!in_expression], collapse = ", "))
}
if (any(!in_annotation & in_expression)) {
  message("  Not found in ATAC gene annotation (skipped): ",
          paste(GENE_PANEL[!in_annotation & in_expression], collapse = ", "))
}
GENE_PANEL <- GENE_PANEL[plottable]
if (length(GENE_PANEL) == 0) stop("No requested genes have both RNA expression and ATAC annotation coordinates.")

# Small symbol-keyed RNA assay for CoveragePlot's expression track -- see
# header for why this is necessary rather than using originalexp directly.
gene_ensembl_ids <- setNames(rownames(fm_rna)[match(GENE_PANEL, fm_rna$gene_symbols)], GENE_PANEL)
rna_panel_data <- GetAssayData(uro, assay = "originalexp", layer = "data")[gene_ensembl_ids, , drop = FALSE]
rownames(rna_panel_data) <- names(gene_ensembl_ids)
uro[["rna_symbols"]] <- CreateAssayObject(data = rna_panel_data)

# ── Shared y-axis scale across all genes ────────────────────────────────────
# CoveragePlot's coverage track auto-scales its y-axis per gene (ymax=NULL
# default), which makes peak heights impossible to compare across genes/
# pages -- each page's axis reflects only its own local max. Signac's
# CoverageTrack() (confirmed by reading its source directly) computes that
# per-gene auto ymax as signif(max(coverage), 2) on a normalized, group-
# scaled, rolling-window-smoothed signal built via CutMatrix() +
# ApplyMatrixByGroup(); reproduced that exact same computation here (via
# Signac's own internal helper functions, not a ggplot/patchwork
# introspection hack -- confirmed those don't reliably expose the rendered
# axis range) for every gene, skipping only the final random downsampling
# step CoverageTrack does purely for plotting-performance -- omitting it
# can only make the computed max equal to or higher than what Signac's own
# downsampled per-gene max would be, never lower, so using it as a single
# shared ymax across all genes is guaranteed not to clip any real peak.
message("\n==> Computing a shared y-axis scale across all genes ...")
compute_gene_covmax <- function(object, gene, assay = "ATAC", group.by = "Age",
                                 window = 100, extend.upstream = 5000, extend.downstream = 5000,
                                 sep = c("-", "-")) {
  region <- Signac:::FindRegion(object = object, region = gene, sep = sep, assay = assay,
                                 extend.upstream = extend.upstream, extend.downstream = extend.downstream)
  cells <- colnames(object)
  cells.per.group <- Signac:::CellsPerGroup(object = object, group.by = group.by)
  obj.groups <- Signac:::GetGroups(object = object, group.by = group.by, idents = NULL)
  obj.groups <- obj.groups[cells]
  reads.per.group <- Signac:::AverageCounts(object = object, assay = assay, group.by = group.by, verbose = FALSE)
  cutmat <- Signac:::CutMatrix(object = object, region = region, assay = assay, cells = cells, verbose = FALSE)
  colnames(cutmat) <- start(region):end(region)
  group.scale.factors <- suppressWarnings(reads.per.group * cells.per.group)
  scale.factor <- median(group.scale.factors)
  coverages <- Signac:::ApplyMatrixByGroup(mat = cutmat, fun = colSums, groups = obj.groups,
                                            group.scale.factors = group.scale.factors,
                                            scale.factor = scale.factor, normalize = TRUE)
  coverages <- group_by(coverages, group)
  coverages <- mutate(coverages, coverage = RcppRoll::roll_sum(x = norm.value, n = window, fill = NA, align = "center"))
  coverages <- ungroup(coverages)
  coverages <- coverages[!is.na(coverages$coverage), ]
  signif(max(coverages$coverage, na.rm = TRUE), digits = 2)
}

Idents(uro) <- "Age"
gene_covmax <- sapply(GENE_PANEL, function(g) {
  v <- tryCatch(compute_gene_covmax(uro, g), error = function(e) NA_real_)
  message(sprintf("  %-8s covmax = %s", g, ifelse(is.na(v), "failed", v)))
  v
})
shared_ymax <- max(gene_covmax, na.rm = TRUE)
message(sprintf("  Shared ymax (max across all genes): %s", shared_ymax))

################################################################################
# Panel G: one CoveragePlot page per gene, common y-axis scale
################################################################################
message("\n==> Building Panel G (per-gene CoveragePlots, all Uro cells by stage) ...")

out_pdf <- file.path(OUT_DIR, "FigG_UrotheliumDevelopment_GenePanel_CoveragePlots.pdf")
pdf(out_pdf, width = 9, height = 8)
for (gene in GENE_PANEL) {
  message("  ", gene)
  p <- tryCatch(
    CoveragePlot(
      object            = uro,
      region            = gene,
      features          = gene,
      assay             = "ATAC",
      expression.assay  = "rna_symbols",
      group.by          = "Age",
      extend.upstream    = 5000,
      extend.downstream  = 5000,
      ymax              = shared_ymax
    ) & patchwork::plot_annotation(
      title = sprintf("%s locus: accessibility and expression by stage", gene)
    ),
    error = function(e) { message("    CoveragePlot failed for ", gene, ": ", e$message); NULL }
  )
  if (!is.null(p)) print(p)
}
dev.off()
message("  Saved: ", out_pdf)

message("\n==> Done.")
