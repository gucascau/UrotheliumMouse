################################################################################
# 06_UrotheliumDevelopment_SCENIC.R
#
# GRN / regulon (SCENIC) analysis of the 715 pseudotime-ordered urothelial
# cells (03_UrotheliumDevelopment_Pseudotime.R), asking:
#   - which TF regulons are active early (E16.5/P0-like pseudotime) and shut
#     off, versus which activate late (mature/barrier-like pseudotime)?
#   - is there a candidate "switch" TF (or small set) whose regulon activity
#     inversely tracks the early/progenitor-associated regulons -- i.e.
#     represses the developmental program while activating the barrier
#     program?
#
# R SCENIC (GENIE3 -> RcisTarget -> AUCell), not pySCENIC -- user's choice.
# mm10 cisTarget databases (mc9nr, refseq-r80, 500bp + 10kb rankings),
# motifs-v9-nr.mgi annotation table, and allTFs_mm.txt were downloaded to
# reference/cisTarget_mm10/ specifically for this analysis (no prior SCENIC
# use in this repo).
#
# Expression input: UrotheliumOnly_rawcounts.mtx/genes.txt/barcodes.txt
# (04a_extract_uro_rawcounts_h5ad.py's true raw-UMI extraction for these 715
# cells -- the RDS objects elsewhere in this pipeline only carry
# log-normalized data, no counts). Gene rownames there are Ensembl IDs;
# mapped to gene_symbols via VisiumLowDevelopmentScripts/output/
# Chen2025_rawcounts/features.tsv.gz, which was built from the SAME source
# h5ad/var table (same 31,671 gene order) for an unrelated Visium
# deconvolution task -- reused here rather than re-deriving the mapping.
#
# Pseudotime: UrotheliumOnly_pseudotime.rds's $Pseudotime column (slingshot,
# rooted fetal -> mature per 03's header comment).
#
# Interpretation anchors: reuses the exact marker panel already established
# in 01_UrotheliumDevelopment_Figures.R's Figure 2 (not redefined here) --
# "Upper-tract developmental" (Pax8, Pax2, Glis3, Fgfr2, Pkhd1, Bicc1) as the
# nephric/ureteric-progenitor axis, "Differentiation / barrier" (Upk1a/b,
# Upk2, Upk3a/b, Krt20, Cldn4, Tjp1) as the mature barrier-program axis.
# Candidate "switch" TFs are regulons whose AUC anti-correlates with the
# progenitor-marker regulons (esp. Pax2/Pax8, if SCENIC recovers a regulon
# for them) and whose own target gene sets are enriched for the barrier
# markers -- ranked, not asserted; SCENIC and the correlation numbers decide.
#
# SCENIC's runSCENIC_* steps write to relative "int/"/"output/" dirs --
# setwd() into a dedicated SCENIC_run/ subdirectory first so these land
# under output/SCENIC_run/{int,output}/ rather than polluting output/
# directly (or worse, wherever the script happened to be launched from).
#
# Input:  UrotheliumDevelopmentScripts/output/UrotheliumOnly_rawcounts.mtx
#         UrotheliumDevelopmentScripts/output/UrotheliumOnly_rawcounts_genes.txt
#         UrotheliumDevelopmentScripts/output/UrotheliumOnly_rawcounts_barcodes.txt
#         UrotheliumDevelopmentScripts/output/UrotheliumOnly_pseudotime.rds
#         VisiumLowDevelopmentScripts/output/Chen2025_rawcounts/features.tsv.gz
#         reference/cisTarget_mm10/*.feather, motifs-v9-nr.mgi*.tbl, allTFs_mm.txt
# Output: UrotheliumDevelopmentScripts/output/SCENIC_run/int/, output/ (SCENIC's
#           own intermediate + regulon/AUC files)
#         UrotheliumOnly_regulonAUC.rds (regulons x cells AUC matrix)
#         UrotheliumOnly_regulon_pseudotime_correlation.csv (per-regulon
#           Spearman rho vs Pseudotime, ranked)
#         Fig9_UrotheliumDevelopment_RegulonActivity_Pseudotime.pdf
#         Fig10_UrotheliumDevelopment_SwitchTF_Candidates.pdf
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(SCENIC)
  library(GENIE3)
  library(RcisTarget)
  library(AUCell)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

SCRIPT_DIR    <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/UrotheliumDevelopmentScripts"
OUT_DIR       <- file.path(SCRIPT_DIR, "output")
SCENIC_RUNDIR <- file.path(OUT_DIR, "SCENIC_run")
dir.create(SCENIC_RUNDIR, showWarnings = FALSE, recursive = TRUE)

DB_DIR       <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/reference/cisTarget_mm10"
FEATURES_MAP <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/VisiumLowDevelopmentScripts/output/Chen2025_rawcounts/features.tsv.gz"

N_CORES <- 8

# ── Marker panel from 01_UrotheliumDevelopment_Figures.R Figure 2 (reused, ──
# ── not redefined) ───────────────────────────────────────────────────────────
DEV_MARKERS     <- c("Pax8", "Pax2", "Glis3", "Fgfr2", "Pkhd1", "Bicc1")
IMMATURE_MARKERS <- c("Trp63", "Krt14", "Mki67", "Top2a")
BARRIER_MARKERS <- c("Upk1a", "Upk1b", "Upk2", "Upk3a", "Upk3b", "Krt20", "Cldn4", "Tjp1")

# ── Load raw counts, map Ensembl -> gene symbol ─────────────────────────────
message("==> Loading raw UMI counts for the 715 Uro-only cells ...")
counts <- Matrix::readMM(file.path(OUT_DIR, "UrotheliumOnly_rawcounts.mtx"))
gene_ids  <- readLines(file.path(OUT_DIR, "UrotheliumOnly_rawcounts_genes.txt"))
barcodes  <- readLines(file.path(OUT_DIR, "UrotheliumOnly_rawcounts_barcodes.txt"))
rownames(counts) <- gene_ids
colnames(counts) <- barcodes
message(sprintf("  %d genes x %d cells", nrow(counts), ncol(counts)))

feat_map <- read.delim(FEATURES_MAP, header = FALSE,
                        col.names = c("gene_id", "gene_symbol", "feature_type"),
                        stringsAsFactors = FALSE)
stopifnot(all(gene_ids %in% feat_map$gene_id))
sym <- feat_map$gene_symbol[match(gene_ids, feat_map$gene_id)]
keep <- !duplicated(sym) & sym != ""
counts <- counts[keep, ]
rownames(counts) <- sym[keep]
message(sprintf("  %d genes after Ensembl->symbol dedup", nrow(counts)))

exprMat <- as.matrix(counts)
rm(counts); gc()

# ── Load Pseudotime metadata, align to the count matrix's cells ────────────
message("==> Loading Pseudotime metadata ...")
pt_obj <- readRDS(file.path(OUT_DIR, "UrotheliumOnly_pseudotime.rds"))
stopifnot(all(colnames(exprMat) %in% colnames(pt_obj)))
cell_meta <- pt_obj@meta.data[colnames(exprMat), c("Age", "Pseudotime")]
rm(pt_obj); gc()

# ── SCENIC setup ─────────────────────────────────────────────────────────────
message("==> Initializing SCENIC ...")
old_wd <- getwd()
setwd(SCENIC_RUNDIR)
on.exit(setwd(old_wd), add = TRUE)

# SCENIC::checkAnnots() (called inside initializeScenic()) hardcodes the
# ranking index column name as "features" regardless of dbIndexCol -- but
# the aertslab-hosted mc9nr feather files ship with that column named
# "motifs" instead, which breaks checkAnnots with a "column doesn't exist"
# error. Fixed at the source: both feather files in DB_DIR were rewritten
# once (interactively, not by this script) with that column renamed
# motifs -> features, so the default dbIndexCol = "features" now works.
#
# RcisTarget also expects a `motifAnnotations_mgi` object to already exist
# in the environment (initializeScenic() -> getDbAnnotations() does
# eval(as.name(paste0("motifAnnotations_", org)))) -- built here from the
# downloaded motifs-v9-nr.mgi tbl rather than relying on a separate
# pre-packaged annotation data package (none installed for mm9/mm10 mgi).
motifAnnotations_mgi <- importAnnotations(file.path(DB_DIR, "motifs-v9-nr.mgi-m0.001-o0.0.tbl"))

dbs <- c(
  "500bp" = "mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather",
  "10kb"  = "mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
)
scenicOptions <- initializeScenic(
  org         = "mgi",
  dbDir       = DB_DIR,
  dbs         = dbs,
  datasetTitle = "UrotheliumDevelopment_Pseudotime",
  nCores      = N_CORES
)
scenicOptions@settings$db_mcVersion <- "v9"
saveRDS(scenicOptions, file.path("int", "scenicOptions.Rds"))

# ── Gene filtering (SCENIC's standard low-expression/low-detection filter) ──
message("==> Filtering genes ...")
genesKept <- geneFiltering(exprMat, scenicOptions = scenicOptions,
                            minCountsPerGene = 3 * 0.01 * ncol(exprMat),
                            minSamples       = ncol(exprMat) * 0.01)
message(sprintf("  %d / %d genes kept", length(genesKept), nrow(exprMat)))
exprMat_filtered <- exprMat[genesKept, ]

runCorrelation(exprMat_filtered, scenicOptions)

# GENIE3 and AUCell both run on log2(counts+1) -- SCENIC's current vignette
# guidance (not the raw counts the very first SCENIC papers used); keeps the
# GRN step and the scoring step on the same expression scale.
exprMat_filtered_log <- log2(exprMat_filtered + 1)

message("==> Running GENIE3 (TF-target GRN inference) ...")
# runGenie3(resumePreviousRun = TRUE) errors out (rather than skipping) once
# its output already exists -- by design, it's telling the caller not to call
# it again, not asking to be called with that flag. A prior run already
# completed GENIE3 (and got as far as building regulons) before failing at
# the AUCell step on a missing doMC package, so skip the call entirely now
# that doMC is installed, since GENIE3's own output didn't change.
genie3_done <- file.exists(getIntName(scenicOptions, "genie3ll"))
if (genie3_done) {
  message("  GENIE3 output already exists (int/1.4_GENIE3_linkList.Rds) -- skipping re-run.")
} else {
  runGenie3(exprMat_filtered_log, scenicOptions, resumePreviousRun = TRUE)
}

# ── Build regulons ───────────────────────────────────────────────────────────
message("==> Building co-expression modules + regulons (RcisTarget) ...")
scenicOptions@settings$nCores <- N_CORES
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethods = c("top5perTarget"))

# ── Score cells (AUCell) ─────────────────────────────────────────────────────
message("==> Scoring cells with AUCell ...")
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered_log)
saveRDS(scenicOptions, file.path("int", "scenicOptions.Rds"))

regulonAUC <- getAUC(loadInt(scenicOptions, "aucell_regulonAUC"))
message(sprintf("  %d regulons scored across %d cells", nrow(regulonAUC), ncol(regulonAUC)))
saveRDS(regulonAUC, file.path(OUT_DIR, "UrotheliumOnly_regulonAUC.rds"))


regulonAUC <- readRDS(file.path(OUT_DIR, "UrotheliumOnly_regulonAUC.rds"))
# scenicOptions <- readRDS(file.path(SCENIC_RUNDIR,"int", "scenicOptions.Rds"))
regulonAUC %>% head()

# regulonTargetsInfo (TF -> target gene table, used below both for regulon
# size and for the barrier-marker target-overlap screen) is also loaded here,
# while cwd is still SCENIC_RUNDIR -- loadInt() resolves scenicOptions'
# relative "int/..." paths against the CURRENT working directory, so this
# must happen before setwd(old_wd) below, not after.
regulonTargets <- loadInt(scenicOptions, "regulonTargetsInfo")
regulon_size <- table(regulonTargets$TF)

setwd(old_wd)

# ── Regulon activity vs pseudotime ───────────────────────────────────────────
message("==> Correlating regulon activity with Pseudotime ...")
regulonAUC <- regulonAUC[, rownames(cell_meta)]  # enforce matching cell order
pseudotime <- cell_meta$Pseudotime
age_levels <- unique(cell_meta$Age)
message(sprintf("  Age strata for concordance check: %s", paste(age_levels, collapse = ", ")))

regulon_cor <- lapply(rownames(regulonAUC), function(reg) {
  auc <- regulonAUC[reg, ]
  valid <- !is.na(pseudotime) & !is.na(auc)
  ct <- suppressWarnings(cor.test(auc[valid], pseudotime[valid], method = "spearman"))
  pooled_rho <- unname(ct$estimate)

  # "Highly conserved within pseudotime" = the pooled trend isn't just an
  # age-batch artifact -- check whether the same sign of correlation holds
  # up within each age stratum on its own (age groups with <10 valid cells
  # are skipped as too noisy to trust a sign from).
  per_age_sign <- sapply(age_levels, function(a) {
    idx <- valid & cell_meta$Age == a
    if (sum(idx) < 10) return(NA_real_)
    ct_age <- suppressWarnings(cor.test(auc[idx], pseudotime[idx], method = "spearman"))
    sign(unname(ct_age$estimate))
  })
  n_age_groups_tested <- sum(!is.na(per_age_sign))
  n_age_groups_concordant <- sum(per_age_sign == sign(pooled_rho), na.rm = TRUE)

  tf_name <- sub("_extended$", "", sub("\\s*\\(.*\\)$", "", reg))
  n_genes <- if (tf_name %in% names(regulon_size)) unname(regulon_size[tf_name]) else NA_integer_
  data.frame(regulon = reg, rho = pooled_rho, p_value = ct$p.value, n_genes = n_genes,
             n_age_groups_concordant = n_age_groups_concordant,
             n_age_groups_tested = n_age_groups_tested)
}) %>% bind_rows()

regulon_cor <- regulon_cor %>%
  mutate(fdr = p.adjust(p_value, method = "BH"),
         frac_age_groups_concordant = n_age_groups_concordant / n_age_groups_tested) %>%
  arrange(rho)

# "TF" alone denotes an activating regulon (extended = TF_extended, more
# permissive motif matches); classify by name for readability in the CSV.
regulon_cor$is_extended <- grepl("_extended", regulon_cor$regulon)
regulon_cor$tf <- sub("_extended$", "", sub("\\s*\\(.*\\)$", "", regulon_cor$regulon))

out_csv <- file.path(OUT_DIR, "UrotheliumOnly_regulon_pseudotime_correlation.csv")
write.csv(regulon_cor, out_csv, row.names = FALSE)
message("  Saved: ", out_csv)
# regulon_cor<-read.csv(out_csv)


message("\n  Top 10 regulons DEcreasing with pseudotime (early-active, shut off):")
print(head(regulon_cor %>% arrange(rho, desc(frac_age_groups_concordant)), 10))
message("\n  Top 10 regulons INcreasing with pseudotime (late-activating):")
print(head(regulon_cor %>% arrange(desc(rho), desc(frac_age_groups_concordant)), 10))

# ── Locate the progenitor-axis reference regulons (Pax2/Pax8 etc, if any) ──
progenitor_regulons <- regulon_cor %>% filter(tf %in% DEV_MARKERS)
message("\n  Progenitor-marker (Fig 2 'Upper-tract developmental') regulons recovered:")
print(progenitor_regulons)

# ── Candidate "switch" TFs: regulons whose AUC anti-correlates with the ────
# ── mean progenitor-regulon AUC AND whose targets enrich for barrier genes ──
message("\n==> Screening for candidate switch TFs ...")
if (nrow(progenitor_regulons) > 0) {
  progenitor_auc <- colMeans(regulonAUC[progenitor_regulons$regulon, , drop = FALSE])

  switch_candidates <- lapply(rownames(regulonAUC), function(reg) {
    auc <- regulonAUC[reg, ]
    ct <- suppressWarnings(cor.test(auc, progenitor_auc, method = "spearman"))
    tf_name <- sub("_extended$", "", sub("\\s*\\(.*\\)$", "", reg))
    targets <- regulonTargets$gene[regulonTargets$TF == tf_name]
    n_barrier_targets <- sum(BARRIER_MARKERS %in% targets)
    data.frame(regulon = reg, tf = tf_name,
               rho_vs_progenitor_regulons = unname(ct$estimate),
               p_value = ct$p.value,
               n_barrier_targets = n_barrier_targets,
               barrier_targets = paste(intersect(BARRIER_MARKERS, targets), collapse = ", "))
  }) %>% bind_rows() %>%
    filter(!tf %in% DEV_MARKERS) %>%
    mutate(fdr = p.adjust(p_value, method = "BH")) %>%
    # frac_age_groups_concordant carried over from regulon_cor as a tiebreaker
    # -- prefer candidates whose own pseudotime trend is conserved across age
    # strata (not just an age-batch artifact) over ones with the same
    # rho/barrier-overlap but a less robust pseudotime trend.
    left_join(regulon_cor %>% select(regulon, frac_age_groups_concordant), by = "regulon") %>%
    arrange(rho_vs_progenitor_regulons, desc(n_barrier_targets), desc(frac_age_groups_concordant))

  out_switch_csv <- file.path(OUT_DIR, "UrotheliumOnly_switchTF_candidates.csv")
  write.csv(switch_candidates, out_switch_csv, row.names = FALSE)
  message("  Saved: ", out_switch_csv)

  message("\n  Top switch-TF candidates (most anti-correlated with progenitor regulons,")
  message("  ranked further by barrier-marker target overlap):")
  print(head(switch_candidates, 15))
} else {
  message("  No regulon recovered for any Fig 2 progenitor marker (Pax8/Pax2/Glis3/Fgfr2/Pkhd1/Bicc1) --")
  message("  cannot anchor the switch-TF screen on them. See UrotheliumOnly_regulon_pseudotime_correlation.csv")
  message("  for the full early-vs-late regulon ranking instead.")
  switch_candidates <- NULL
}

# ── Plots ─────────────────────────────────────────────────────────────────────
message("\n==> Plotting ...")

# Fig 9: heatmap-style view -- top 15 early + top 15 late regulons' AUC vs
# pseudotime (binned into deciles), analogous in spirit to Fig6's smoothed
# marker-expression-vs-pseudotime plot but for regulon activity.
top_early <- head(regulon_cor %>% filter (grepl("_extended",regulon)) %>% arrange(rho, desc(frac_age_groups_concordant)), 10)$regulon
top_late  <- head(regulon_cor %>% filter (grepl("_extended",regulon))%>% arrange(desc(rho), desc(frac_age_groups_concordant)), 10)$regulon
plot_regulons <- unique(c(top_early, top_late))

pt_bins <- cut(pseudotime, breaks = 10, labels = FALSE)
auc_long <- as.data.frame(t(regulonAUC[plot_regulons, , drop = FALSE])) %>%
  mutate(cell = rownames(.), pseudotime_bin = pt_bins) %>%
  pivot_longer(-c(cell, pseudotime_bin), names_to = "regulon", values_to = "AUC") %>%
  group_by(regulon, pseudotime_bin) %>%
  summarise(mean_AUC = mean(AUC, na.rm = TRUE), .groups = "drop") %>%
  mutate(direction = ifelse(regulon %in% top_early, "Early -> shuts off", "Late-activating"))

auc_long %>% arrange(desc(mean_AUC))
fig9 <- ggplot(auc_long, aes(x = pseudotime_bin, y = reorder(regulon, log2(mean_AUC+0.01)), fill = log2(mean_AUC+0.01))) +
  geom_tile() +
  facet_grid(direction ~ ., scales = "free_y", space = "free_y") +
  scale_fill_viridis_c(option = "plasma", name = "Log2 mean\nregulon AUC") +
  labs(x = "Pseudotime (decile bin)", y = NULL,
       title = "Figure 9. Regulon activity across pseudotime",
       subtitle = "Top 10 early-shutoff and top 10 late-activating regulons") +
  theme_minimal(base_size = 10) +
  theme(strip.text.y = element_text(angle = 0))
ggsave(file.path(OUT_DIR, "Fig9_UrotheliumDevelopment_RegulonActivity_Pseudotime.pdf"),
       plot = fig9, width = 7, height = 6)
message("  Saved: Fig9_UrotheliumDevelopment_RegulonActivity_Pseudotime.pdf")

# Fig 10: switch-TF candidate regulons vs the progenitor-regulon signature,
# scatter over pseudotime.
if (!is.null(switch_candidates) && nrow(switch_candidates) > 0) {
  top_switch <- head(switch_candidates$regulon, 4)
  switch_df <- data.frame(
    pseudotime = pseudotime,
    progenitor_AUC = progenitor_auc
  )
  for (reg in top_switch) switch_df[[reg]] <- regulonAUC[reg, ]
  switch_long <- switch_df %>%
    pivot_longer(-pseudotime, names_to = "regulon", values_to = "AUC")

  fig10 <- ggplot(switch_long, aes(x = pseudotime, y = AUC, color = regulon)) +
    geom_point(alpha = 0.4, size = 0.8) +
    geom_smooth(se = FALSE, method = "loess") +
    labs(x = "Pseudotime", y = "Regulon AUC",
         title = "Figure 10. Candidate switch-TF regulons vs. progenitor-regulon signature",
         subtitle = paste("Progenitor signature = mean AUC of", paste(progenitor_regulons$tf, collapse = ", "), "regulons")) +
    theme_minimal(base_size = 10)
  ggsave(file.path(OUT_DIR, "Fig10_UrotheliumDevelopment_SwitchTF_Candidates.pdf"),
         plot = fig10, width = 9, height = 6)
  message("  Saved: Fig10_UrotheliumDevelopment_SwitchTF_Candidates.pdf")
}

message("\n==> Done.")
