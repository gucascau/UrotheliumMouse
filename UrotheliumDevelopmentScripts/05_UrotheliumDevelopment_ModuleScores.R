################################################################################
# 05_UrotheliumDevelopment_ModuleScores.R
#
# Collapse the individual-gene marker panel from Figure 2/6 (plus two new
# gene sets) into 5 interpretable per-cell module scores via Seurat's
# AddModuleScore (expression averaged per cell, minus a matched-expression-
# bin control set -- standard control-corrected module scoring, not a raw
# mean). Scored on the full 31,671-gene "originalexp" assay (Ensembl-ID
# rownames) so AddModuleScore has enough background genes to bin controls
# from; gene symbols below are mapped to Ensembl IDs via meta.features, same
# pattern as 03/04.
#
# Module definitions:
#   Upper-tract identity            -- Pax8, Pax2, Glis3, Fgfr2, Pkhd1, Bicc1
#                                       (same panel as Figure 2/6; Fig 2's
#                                       gene-by-gene recheck showed this is
#                                       NOT uniformly "fetal" -- Pkhd1/Glis3
#                                       actually rise with age -- so this
#                                       module's net trend is a genuine
#                                       question, not an assumed decline)
#   Barrier / differentiation       -- Upk1a/1b/2/3a/3b, Krt20, Cldn4, Tjp1
#   Basement membrane / epithelial
#     organization                  -- Col4a3/4/5 (type IV collagen), Lama5,
#                                       Nid1 (basement membrane structure),
#                                       Magi1, Dsp (junctional organization)
#   Proliferation / immaturity      -- Mki67, Top2a, Pcna, Ccnb1, Cdk1 (pure
#                                       cell-cycle genes only -- Trp63/Krt14
#                                       deliberately excluded: the Figure 2
#                                       recheck showed these track a
#                                       persistent basal population that
#                                       rises postnatally and is NOT an
#                                       immaturity signal)
#   Stress / remodeling             -- Fos, Jun, Junb (immediate-early/
#                                       stress response), Spp1, Fn1, Timp1
#                                       (ECM remodeling) -- included per
#                                       request as exploratory ("if present");
#                                       report the result honestly even if it
#                                       doesn't show a clean stage trend
#
# Figure 8. Module score dynamics
#   Violin + boxplot of each module score by developmental stage, with the
#   Spearman correlation (module score vs. chronological stage) annotated
#   per panel.
#
# Input:  UrotheliumDevelopmentScripts/output/UrotheliumOnly_pseudotime_cytotrace.rds
#         (written by 04_UrotheliumDevelopment_CytoTRACE.R)
# Output: Fig8_UrotheliumDevelopment_ModuleScores.pdf
#         UrotheliumOnly_full_annotated.rds
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(viridisLite)
})

SCRIPT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/UrotheliumDevelopmentScripts"
OUT_DIR    <- file.path(SCRIPT_DIR, "output")

STAGE_ORDER  <- c("E16.5", "P0", "W3", "W12", "W52", "W92")
STAGE_COLORS <- setNames(viridisLite::viridis(length(STAGE_ORDER)), STAGE_ORDER)

# ── Load ─────────────────────────────────────────────────────────────────────
message("==> Loading Uro-only object ...")
in_rds <- file.path(OUT_DIR, "UrotheliumOnly_pseudotime_cytotrace.rds")
if (!file.exists(in_rds)) {
  stop("Missing ", in_rds, " -- run 04_UrotheliumDevelopment_CytoTRACE.R first.")
}
uro <- readRDS(in_rds)
uro$Age <- factor(uro$Age, levels = STAGE_ORDER)
DefaultAssay(uro) <- "originalexp"

################################################################################
# Module scores
################################################################################
MODULE_ORDER <- c("Upper-tract identity", "Proliferation / immaturity",
                   "Basement membrane / epithelial organization",
                   "Barrier / differentiation", "Stress / remodeling")

ModuleSets_symbols <- list(
  "Upper-tract identity" = c("Pax8", "Pax2", "Glis3", "Fgfr2", "Pkhd1", "Bicc1"),
  "Proliferation / immaturity" = c("Mki67", "Top2a", "Pcna", "Ccnb1", "Cdk1"),
  "Basement membrane / epithelial organization" =
    c("Col4a3", "Col4a4", "Col4a5", "Lama5", "Nid1", "Magi1", "Dsp"),
  "Barrier / differentiation" =
    c("Upk1a", "Upk1b", "Upk2", "Upk3a", "Upk3b", "Krt20", "Cldn4", "Tjp1"),
  "Stress / remodeling" = c("Fos", "Jun", "Junb", "Spp1", "Fn1", "Timp1")
)
ModuleSets_symbols <- ModuleSets_symbols[MODULE_ORDER]

fm <- uro[["originalexp"]]@meta.features
ModuleSets_ens <- lapply(ModuleSets_symbols, function(syms) {
  idx <- match(syms, fm$gene_symbols)
  if (any(is.na(idx))) {
    stop("Missing gene symbol(s): ", paste(syms[is.na(idx)], collapse = ", "))
  }
  rownames(fm)[idx]
})

message("\n==> Computing module scores (Seurat::AddModuleScore) ...")
uro <- AddModuleScore(uro, features = ModuleSets_ens, name = "Module",
                       assay = "originalexp", seed = 1)

# AddModuleScore names columns Module1..Module5 in list order -- rename to
# the actual module names.
module_cols <- paste0("Module", seq_along(MODULE_ORDER))
stopifnot(all(module_cols %in% colnames(uro@meta.data)))
colnames(uro@meta.data)[match(module_cols, colnames(uro@meta.data))] <- MODULE_ORDER

message("\n  Mean module score by stage:")
print(uro@meta.data %>% group_by(Age) %>%
        summarise(across(all_of(MODULE_ORDER), mean), .groups = "drop"))

stage_rank <- as.integer(uro$Age)
rho_module <- sapply(MODULE_ORDER, function(m) {
  cor(uro@meta.data[[m]], stage_rank, method = "spearman")
})
message("\n  Spearman correlation(module score, chronological stage):")
print(round(rho_module, 3))

saveRDS(uro, file.path(OUT_DIR, "UrotheliumOnly_full_annotated.rds"))
message("\n  Saved: UrotheliumOnly_full_annotated.rds")

################################################################################
# Figure 8: module score dynamics
################################################################################
message("\n==> Building Figure 8 (module score dynamics) ...")

module_long <- uro@meta.data %>%
  select(Age, all_of(MODULE_ORDER)) %>%
  pivot_longer(-Age, names_to = "Module", values_to = "Score") %>%
  mutate(Module = factor(Module, levels = MODULE_ORDER))

rho_labels <- data.frame(
  Module = factor(MODULE_ORDER, levels = MODULE_ORDER),
  label = sprintf("rho = %.2f", rho_module)
)

fig8 <- ggplot(module_long, aes(x = Age, y = Score, fill = Age)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  geom_text(data = rho_labels, aes(x = -Inf, y = Inf, label = label),
            inherit.aes = FALSE, hjust = -0.1, vjust = 1.4, size = 3.3) +
  scale_fill_manual(values = STAGE_COLORS, guide = "none") +
  facet_wrap(~ Module, ncol = 3, scales = "free_y") +
  labs(x = "Developmental stage", y = "Module score",
       title = "Figure 8. Module score dynamics across development") +
  theme_classic(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 9),
        plot.title = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUT_DIR, "Fig8_UrotheliumDevelopment_ModuleScores.pdf"),
       fig8, width = 12, height = 8)
message("  Saved: Fig8_UrotheliumDevelopment_ModuleScores.pdf")

message("\n==> Done.")
