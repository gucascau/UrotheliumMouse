################################################################################
# 01_UrotheliumDevelopment_Figures.R
#
# Developmental atlas figures for the Chen 2025 NatGenet sex-specific lifespan
# mouse kidney multiome object (whole-kidney, E16.5 -> W92, snRNA).
#
# Figure 1. Developmental atlas overview
#   Panel A: integrated UMAP colored by developmental stage (Age).
#   Panel B: integrated UMAP colored by a 4-class lineage grouping that
#            foregrounds the urothelium and its developmentally-related
#            (ureteric bud-derived) epithelium against nephron epithelium and
#            non-epithelial (stroma/vascular/immune) cells.
#
# Figure 2. Marker dynamics across development (Uro cells only)
#   Dot plot of upper-tract developmental / immature-proliferative /
#   differentiation-barrier marker sets, grouped by developmental stage, to
#   show the fetal-like -> mature epithelial identity shift.
#
# Figure 3. Urothelium proportion across development
#   Per-sample % of cells classified as Uro at each stage (bar = mean +/- SEM,
#   points = individual samples), showing how the urothelial fraction of the
#   captured kidney changes over the lifespan series.
#
# Notes on the input object:
#   - It is a Seurat object (despite the "_zellkonvertedConverted" filename),
#     single assay "originalexp", single layer "data" (already log-normalized
#     -- no counts layer, so do NOT re-run NormalizeData).
#   - rownames(obj) are Ensembl IDs; gene symbols live in
#     obj[["originalexp"]]@meta.features$gene_symbols. We build a small
#     symbol-keyed "Markers" assay for the marker panel used in Figure 2
#     rather than renaming all ~31,700 features.
#   - Precomputed reductions X_pca / X_umap are used as-is (no re-clustering).
#   - "Age" is the developmental-stage column: E16.5, P0, W3, W12, W52, W92.
#   - "celltype_final" has 26 levels; the urothelial population is "Uro"
#     (n = 75/94/127/159/125/135 cells at E16.5/P0/W3/W12/W52/W92).
#
# Input:  MultiOmicSpatialMouseKidney_SexSpecific_lifespan_Chen2025NatGenet_zellkonvertedConverted.rds
# Output: Fig1_UrotheliumDevelopment_UMAP_Overview.pdf
#         Fig2_UrotheliumDevelopment_MarkerDynamics_DotPlot.pdf
#         Fig3_UrotheliumDevelopment_ProportionBySample.pdf
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(viridisLite)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
RAW_BASE   <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/RawMouseSingleCellDatasets"
SCRIPT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/UrotheliumDevelopmentScripts"
OUT_DIR    <- file.path(SCRIPT_DIR, "output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OBJECT_RDS <- file.path(RAW_BASE,
  "MultiOmicSpatialMouseKidney_SexSpecific_lifespan_Chen2025NatGenet_zellkonvertedConverted.rds")

# ── Load object ─────────────────────────────────────────────────────────────
message("==> Loading Chen2025 developmental atlas object ...")
obj <- readRDS(OBJECT_RDS)
message(sprintf("  %d cells x %d genes", ncol(obj), nrow(obj)))
DefaultAssay(obj) <- "originalexp"

# ── Developmental stage ordering ────────────────────────────────────────────
STAGE_ORDER <- c("E16.5", "P0", "W3", "W12", "W52", "W92")
stopifnot(setequal(unique(obj$Age), STAGE_ORDER))
obj$Age <- factor(obj$Age, levels = STAGE_ORDER)
STAGE_COLORS <- setNames(viridisLite::viridis(length(STAGE_ORDER)), STAGE_ORDER)

# ── 4-class lineage grouping for Figure 1 Panel B ───────────────────────────
# Grouping reflects kidney developmental biology: the urothelium and the
# collecting-duct system both derive from the ureteric bud lineage, whereas
# nephron epithelium (PT/LOH/DCT/PEC/Podo + their progenitors) arises
# separately from nephron progenitors via mesenchymal-to-epithelial
# transition. Everything else (stroma/vascular/immune) is collapsed so the
# legend stays readable and the urothelial/epithelial signal is the focus.
NEPHRON_EPI <- c("NP", "NP_proliferate", "PT", "PT(S1)", "PT(S2)", "PT(S3)",
                 "LOH_AL", "LOH_AL_proliferating", "LOH_DL", "DCT", "PEC", "Podo")
COLLECTING_EPI <- c("UBP", "CNT", "CD_PC", "CD_PC_Mix", "CD_IC", "IM", "IM_proliferate")

CELLCLASS_LEVELS <- c(
  "Urothelium",
  "Collecting system / ureteric lineage",
  "Nephron epithelium",
  "Non-epithelial (stroma / vascular / immune)"
)

obj$CellClass <- dplyr::case_when(
  obj$celltype_final == "Uro"            ~ CELLCLASS_LEVELS[1],
  obj$celltype_final %in% COLLECTING_EPI ~ CELLCLASS_LEVELS[2],
  obj$celltype_final %in% NEPHRON_EPI    ~ CELLCLASS_LEVELS[3],
  TRUE                                   ~ CELLCLASS_LEVELS[4]
)
obj$CellClass <- factor(obj$CellClass, levels = CELLCLASS_LEVELS)
message("CellClass composition:")
print(table(obj$CellClass))

# Okabe-Ito colorblind-safe palette; urothelium in the warm/bold color so it
# pops against the two related-but-distinct epithelial classes and the grey
# non-epithelial background.
CELLCLASS_COLORS <- setNames(
  c("#D55E00", "#0072B2", "#009E73", "#B3B3B3"),
  CELLCLASS_LEVELS
)

################################################################################
# Figure 1: Developmental atlas overview
################################################################################
message("\n==> Building Figure 1 (UMAP overview) ...")

## raster=FALSE: Seurat's default >100k-cell rasterization (scattermore)
## alpha-blends heavily overlapping points in dense UMAP clusters, washing
## the fill out to near-white. Plotting as vectors keeps colors solid; the
## resulting PDF (~10-20 MB) is still a normal file size for 203k cells.
p1a <- DimPlot(obj, reduction = "X_umap", group.by = "Age",
               cols = STAGE_COLORS, pt.size = 0.3, raster = FALSE) +
  ggtitle("Developmental stage") +
  labs(color = "Stage") +
  coord_fixed() +
  theme(aspect.ratio = 1, plot.title = element_text(face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 4)))

## Note: Seurat's DimPlot `order=` argument requires the CRAN package
## ggrastr (not installed here) once the default >100k-cell rasterization
## kicks in. Instead we pre-sort cells within the object so Urothelium rows
## are drawn last (i.e. on top) without needing ordered rasterization.
draw_rank <- match(obj$CellClass, CELLCLASS_LEVELS)  # Urothelium = 1 ... Non-epithelial = 4
draw_order <- order(draw_rank, decreasing = TRUE)     # Non-epithelial first, Urothelium last (on top)
obj_p1b <- obj[, draw_order]

p1b <- DimPlot(obj_p1b, reduction = "X_umap", group.by = "CellClass",
               cols = CELLCLASS_COLORS,
               pt.size = 0.3, raster = FALSE) +
  ggtitle("Urothelial & epithelial lineages") +
  labs(color = NULL) +
  coord_fixed() +
  theme(aspect.ratio = 1, plot.title = element_text(face = "bold"),
        legend.text = element_text(size = 8)) +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 1))

fig1 <- (p1a | p1b) +
  plot_annotation(title = "Figure 1. Developmental atlas overview",
                   theme = theme(plot.title = element_text(face = "bold", size = 14)))

ggsave(file.path(OUT_DIR, "Fig1_UrotheliumDevelopment_UMAP_Overview.pdf"),
       fig1, width = 12, height = 5.5)
message("  Saved: Fig1_UrotheliumDevelopment_UMAP_Overview.pdf")

################################################################################
# Figure 2: Marker dynamics across development (Uro cells only)
################################################################################
message("\n==> Building Figure 2 (marker dynamics dot plot) ...")

uro <- subset(obj, subset = celltype_final == "Uro")
uro <- droplevels(uro)
message("  Uro cells per stage:")
print(table(uro$Age))

# Save the clean subset (before the marker-panel "Markers" assay below is
# attached) for reuse by 02_UrotheliumDevelopment_Subcluster_UMAP.R.
saveRDS(uro, file.path(OUT_DIR, "UrotheliumOnly_raw.rds"))
message("  Saved: UrotheliumOnly_raw.rds")

# ── Marker panel (symbols); left-to-right narrative = fetal -> mature ──────
FeatureSets <- list(
  "Upper-tract developmental" = c("Pax8", "Pax2", "Glis3", "Fgfr2", "Pkhd1", "Bicc1"),
  "Immature / proliferative"  = c("Trp63", "Krt14", "Mki67", "Top2a"),
  "Differentiation / barrier" = c("Upk1a", "Upk1b", "Upk2", "Upk3a", "Upk3b",
                                   "Krt20", "Cldn4", "Tjp1")
)
marker_symbols <- unlist(FeatureSets, use.names = FALSE)

# Map symbols -> Ensembl IDs via meta.features, then build a small
# symbol-keyed assay (rownames(obj) are Ensembl IDs in the full object).
fm  <- obj[["originalexp"]]@meta.features
idx <- match(marker_symbols, fm$gene_symbols)
if (any(is.na(idx))) {
  stop("Marker symbols not found in gene_symbols: ",
       paste(marker_symbols[is.na(idx)], collapse = ", "))
}
ens_ids <- rownames(fm)[idx]

mat <- GetAssayData(uro, assay = "originalexp", layer = "data")[ens_ids, , drop = FALSE]
rownames(mat) <- marker_symbols
uro[["Markers"]] <- CreateAssayObject(data = mat)
DefaultAssay(uro) <- "Markers"

fig2 <- DotPlot(uro, features = FeatureSets, group.by = "Age", assay = "Markers",
                cols = c("lightgrey", "red"), dot.scale = 8) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text  = element_text(face = "bold", size = 9),
        plot.title  = element_text(face = "bold", size = 14)) +
  guides(
    color = guide_colorbar(title = "AveExp"),
    size  = guide_legend(title = "PerExp")
  ) +
  labs(x = NULL, y = "Developmental stage",
       title = "Figure 2. Marker dynamics across development (Uro cells)")

ggsave(file.path(OUT_DIR, "Fig2_UrotheliumDevelopment_MarkerDynamics_DotPlot.pdf"),
       fig2, width = 10, height = 4.5)
message("  Saved: Fig2_UrotheliumDevelopment_MarkerDynamics_DotPlot.pdf")

################################################################################
# Figure 3: Urothelium proportion across development
################################################################################
message("\n==> Building Figure 3 (urothelium proportion by stage) ...")

# Proportion per sample_id (not pooled across the whole stage) so the plot
# shows the actual biological/technical spread at each stage, not just one
# number -- each sample_id maps to exactly one Age (checked below).
sample_age_map <- obj@meta.data %>% distinct(sample_id, Age)
stopifnot(!any(duplicated(sample_age_map$sample_id)))

prop_df <- obj@meta.data %>%
  group_by(sample_id, Age) %>%
  summarise(n_total = n(),
            n_uro   = sum(celltype_final == "Uro"),
            .groups = "drop") %>%
  mutate(pct_uro = 100 * n_uro / n_total)

message("  Uro proportion per sample (%):")
print(as.data.frame(prop_df %>% arrange(Age)))

stage_summary <- prop_df %>%
  group_by(Age) %>%
  summarise(mean_pct  = mean(pct_uro),
            sem_pct   = sd(pct_uro) / sqrt(n()),
            n_samples = n(),
            .groups = "drop")

fig3 <- ggplot(stage_summary, aes(x = Age, y = mean_pct, fill = Age)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean_pct - sem_pct, ymax = mean_pct + sem_pct),
                width = 0.15, linewidth = 0.4) +
  geom_jitter(data = prop_df, aes(x = Age, y = pct_uro),
              inherit.aes = FALSE, width = 0.1, size = 1.8, alpha = 0.7) +
  scale_fill_manual(values = STAGE_COLORS, guide = "none") +
  labs(x = "Developmental stage", y = "Urothelium (% of cells per sample)",
       title = "Figure 3. Urothelium proportion across development",
       subtitle = "Bars = mean +/- SEM across samples; points = individual samples") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14))

ggsave(file.path(OUT_DIR, "Fig3_UrotheliumDevelopment_ProportionBySample.pdf"),
       fig3, width = 6, height = 4.5)
message("  Saved: Fig3_UrotheliumDevelopment_ProportionBySample.pdf")

message("\n==> Done.")
