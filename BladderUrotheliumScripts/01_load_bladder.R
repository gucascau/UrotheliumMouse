################################################################################
# Script 01: Load Bladder Urothelium Datasets → individual Seurat RDS
#
# Datasets:
# ┌──────────────────────┬─────────────┬────────────┬──────────────────────────┐
# │ sample_id            │ Accession   │ Technology │ Condition                │
# ├──────────────────────┼─────────────┼────────────┼──────────────────────────┤
# │ BladderNormal1       │ GSM3723360  │ 10X_v2    │ Normal bladder           │
# │ BladderNormal2       │ GSM3723361  │ 10X_v2    │ Normal bladder           │
# │ BladderUro1          │ GSM4970435  │ 10X       │ Healthy urothelium       │
# │ BladderWT1           │ GSM4201633  │ 10X       │ WT mouse bladder         │
# │ BladderB8W1          │ GSM5014059  │ 10X       │ Healthy 8 wks            │
# │ BladderB8W2          │ GSM5014060  │ 10X       │ Healthy 8 wks            │
# │ BladderH48_1         │ GSM5014061  │ 10X       │ Acute CPP 48h            │
# │ BladderH48_2         │ GSM5014062  │ 10X       │ Acute CPP 48h            │
# │ BladderD11_1         │ GSM5014063  │ 10X       │ Chronic CPP 11 days      │
# │ BladderD11_2         │ GSM5014064  │ 10X       │ Chronic CPP 11 days      │
# └──────────────────────┴─────────────┴────────────┴──────────────────────────┘
#
# Format notes:
#   GSE129845_mouse  – standard 10X MTX (genes.tsv.gz, gene_col=2)
#   GSE163029        – standard 10X MTX (genes.tsv.gz, gene_col=2)
#   GSM4201633       – 10X v3 MTX (features.tsv.gz, 3 cols); contains mixed
#                      mouse+human features — filter to ENSMUSG only
#   GSE164557        – dense ASCII tab-delimited XLS (Gene_ID | Symbol | cells)
################################################################################

library(Matrix)
library(data.table)
library(Seurat)

DATA_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/BladderUrothelium"
RAW_DIR  <- file.path(DATA_DIR, "Raw")
OUT_DIR  <- file.path(DATA_DIR, "seurat_objects")
dir.create(OUT_DIR, showWarnings = FALSE)

################################################################################
# Master metadata table
################################################################################

sample_metadata <- data.frame(
  sample_id  = c("BladderNormal1", "BladderNormal2",
                 "BladderUro1",
                 "BladderWT1",
                 "BladderB8W1",    "BladderB8W2",
                 "BladderH48_1",   "BladderH48_2",
                 "BladderD11_1",   "BladderD11_2"),
  gsm_id     = c("GSM3723360", "GSM3723361",
                 "GSM4970435",
                 "GSM4201633",
                 "GSM5014059", "GSM5014060",
                 "GSM5014061", "GSM5014062",
                 "GSM5014063", "GSM5014064"),
  technology = c("10X_v2", "10X_v2",
                 "10X",
                 "10X",
                 "10X", "10X",
                 "10X", "10X",
                 "10X", "10X"),
  species    = rep("mouse", 10),
  tissue     = rep("bladder", 10),
  condition  = c("Normal", "Normal",
                 "Healthy_Urothelium",
                 "WT_Bladder",
                 "Healthy_8wk",    "Healthy_8wk",
                 "Acute_CPP_48h",  "Acute_CPP_48h",
                 "Chronic_CPP_11d","Chronic_CPP_11d"),
  paper      = c("PMID_31462402", "PMID_31462402",
                 "PMID_33538002",
                 "GSE141348",
                 "GSE164557", "GSE164557",
                 "GSE164557", "GSE164557",
                 "GSE164557", "GSE164557"),
  has_annotation = rep(FALSE, 10),
  stringsAsFactors = FALSE
)

get_meta <- function(sid) sample_metadata[sample_metadata$sample_id == sid, ]


################################################################################
# Helpers
################################################################################

load_10x_mtx <- function(barcode_file, gene_file, matrix_file, gene_col = 2) {
  bc_con <- file(barcode_file, "rb")
  barcodes <- readLines(gzcon(bc_con))
  close(bc_con)

  genes_df   <- read.table(gzfile(gene_file), header = FALSE, sep = "\t",
                            stringsAsFactors = FALSE)
  gene_names <- make.unique(genes_df[, gene_col])

  mat_con <- file(matrix_file, "rb")
  mat     <- readMM(gzcon(mat_con))
  close(mat_con)

  rownames(mat) <- gene_names
  colnames(mat) <- barcodes
  mat
}

make_seurat <- function(mat, sample_id, meta_row,
                        min_cells = 3, min_features = 200) {
  so <- CreateSeuratObject(counts = mat, project = sample_id,
                           min.cells = min_cells, min.features = min_features)
  so$sample_id      <- sample_id
  so$gsm_id         <- meta_row$gsm_id
  so$technology     <- meta_row$technology
  so$species        <- meta_row$species
  so$tissue         <- meta_row$tissue
  so$condition      <- meta_row$condition
  so$paper          <- meta_row$paper
  so$has_annotation <- meta_row$has_annotation
  so
}


################################################################################
# SECTION 1: GSM3723360 + GSM3723361 — GSE129845_mouse, 10X v2
# Yu et al. 2019 (PMID_31462402): human & mouse bladder single-cell map
################################################################################

gse129845_samples <- list(
  BladderNormal1 = list(
    b = "GSM3723360_Sample4_barcodes.tsv.gz",
    g = "GSM3723360_Sample4_genes.tsv.gz",
    m = "GSM3723360_Sample4_matrix.mtx.gz"
  ),
  BladderNormal2 = list(
    b = "GSM3723361_Sample5_barcodes.tsv.gz",
    g = "GSM3723361_Sample5_genes.tsv.gz",
    m = "GSM3723361_Sample5_matrix.mtx.gz"
  )
)

for (sid in names(gse129845_samples)) {
  f   <- gse129845_samples[[sid]]
  dir <- file.path(RAW_DIR, "GSE129845_mouse")
  message(sprintf("Loading %s (%s)...", sid, get_meta(sid)$gsm_id))
  mat <- load_10x_mtx(file.path(dir, f$b), file.path(dir, f$g),
                      file.path(dir, f$m))
  so  <- make_seurat(mat, sid, get_meta(sid))
  saveRDS(so, file.path(OUT_DIR, paste0(sid, "_seurat.rds")))
  message(sprintf("  %d cells × %d genes  →  saved", ncol(so), nrow(so)))
  rm(mat, so); gc()
}


################################################################################
# SECTION 2: GSM4970435 — GSE163029, 10X
# Liu et al. 2021 (PMID_33538002): 8 urothelial subpopulations
################################################################################

message("Loading BladderUro1 (GSM4970435)...")
mat <- load_10x_mtx(
  file.path(RAW_DIR, "GSE163029", "GSM4970435_barcodes.tsv.gz"),
  file.path(RAW_DIR, "GSE163029", "GSM4970435_genes.tsv.gz"),
  file.path(RAW_DIR, "GSE163029", "GSM4970435_matrix.mtx.gz")
)
so  <- make_seurat(mat, "BladderUro1", get_meta("BladderUro1"))
saveRDS(so, file.path(OUT_DIR, "BladderUro1_seurat.rds"))
message(sprintf("  %d cells × %d genes  →  saved", ncol(so), nrow(so)))
rm(mat, so); gc()


################################################################################
# SECTION 3: GSM4201633 — 10X v3 features file, barnyard (mouse + human)
# This sample was aligned to a combined mm10 + GRCh38 reference.
# Mouse features: col1 = "mm10___ENSMUSG...", col2 = "mm10___Xkr4"
# Human features: col1 = "GRCh38_ENSG...",   col2 = "GRCh38_MIR1302-2HG"
# Strategy: keep only mm10___ rows, then strip the "mm10___" prefix from
# col2 so the resulting gene symbols ("Xkr4", "Sox17", ...) match the other
# samples.
################################################################################

message("Loading BladderWT1 (GSM4201633)...")
feat_path <- file.path(RAW_DIR, "GSM4201633",
                       "GSM4201633_wt_mouse_bladder_features.tsv.gz")
feat_df   <- read.table(gzfile(feat_path), header = FALSE, sep = "\t",
                         stringsAsFactors = FALSE)
# Col 1: "mm10___ENSMUSG..." or "GRCh38_ENSG..."
# Col 2: "mm10___Xkr4"      or "GRCh38_MIR1302-2HG"
# Col 3: feature type (all "Gene Expression" here)
mouse_ge_idx <- which(grepl("^mm10___", feat_df[, 1]) &
                      feat_df[, 3] == "Gene Expression")
message(sprintf("  Features total: %d  |  mm10 Gene Expression: %d  |  GRCh38: %d",
                nrow(feat_df), length(mouse_ge_idx),
                sum(grepl("^GRCh38_", feat_df[, 1]))))

bc_con   <- file(file.path(RAW_DIR, "GSM4201633",
                            "GSM4201633_wt_mouse_bladder_barcodes.tsv.gz"), "rb")
barcodes <- readLines(gzcon(bc_con))
close(bc_con)

mat_con  <- file(file.path(RAW_DIR, "GSM4201633",
                            "GSM4201633_wt_mouse_bladder_matrix.mtx.gz"), "rb")
mat_full <- readMM(gzcon(mat_con))
close(mat_con)
colnames(mat_full) <- barcodes

mat_mouse           <- mat_full[mouse_ge_idx, ]
# Strip "mm10___" prefix so names match all other datasets (e.g. "Xkr4")
clean_symbols       <- sub("^mm10___", "", feat_df[mouse_ge_idx, 2])
rownames(mat_mouse) <- make.unique(clean_symbols)
rm(feat_df, mat_full, barcodes, clean_symbols); gc()

so <- make_seurat(mat_mouse, "BladderWT1", get_meta("BladderWT1"))
saveRDS(so, file.path(OUT_DIR, "BladderWT1_seurat.rds"))
message(sprintf("  %d cells × %d genes  →  saved", ncol(so), nrow(so)))
rm(mat_mouse, so); gc()


################################################################################
# SECTION 4: GSE164557 — 6 samples, dense ASCII tab-delimited (.xls extension)
# File format: Gene_ID \t Symbol \t cell1 \t cell2 \t ...
# Conditions: B8W = healthy 8 wks; H48 = acute CPP 48h; D11 = chronic CPP 11d
################################################################################

xls_samples <- list(
  BladderB8W1  = "GSM5014059_B8W-1_gene_cell_exprs_table.xls",
  BladderB8W2  = "GSM5014060_B8W-2_gene_cell_exprs_table.xls",
  BladderH48_1 = "GSM5014061_H48-1_gene_cell_exprs_table.xls",
  BladderH48_2 = "GSM5014062_H48-2_gene_cell_exprs_table.xls",
  BladderD11_1 = "GSM5014063_D11-1_gene_cell_exprs_table.xls",
  BladderD11_2 = "GSM5014064_D11-2_gene_cell_exprs_table.xls"
)

for (sid in names(xls_samples)) {
  fname <- xls_samples[[sid]]
  fpath <- file.path(RAW_DIR, "GSE164557", fname)
  message(sprintf("Loading %s (%s)...", sid, fname))

  # Dense tab-delimited: col1 = Ensembl ID, col2 = Symbol, rest = cells
  dt        <- fread(fpath, header = TRUE, sep = "\t", data.table = FALSE,
                     check.names = FALSE)
  gene_syms <- make.unique(dt[[2]])
  cells     <- colnames(dt)[-(1:2)]
  mat_dense <- as.matrix(dt[, -(1:2), drop = FALSE])
  rownames(mat_dense) <- gene_syms
  colnames(mat_dense) <- cells
  mat_sparse <- as(mat_dense, "dgCMatrix")
  rm(dt, mat_dense); gc()

  so <- make_seurat(mat_sparse, sid, get_meta(sid))
  saveRDS(so, file.path(OUT_DIR, paste0(sid, "_seurat.rds")))
  message(sprintf("  %d cells × %d genes  →  saved", ncol(so), nrow(so)))
  rm(mat_sparse, so); gc()
}


################################################################################
# Verification summary
################################################################################

rds_files  <- sort(list.files(OUT_DIR, pattern = "_seurat\\.rds$",
                               full.names = TRUE))
verify_df  <- do.call(rbind, lapply(rds_files, function(f) {
  sid <- sub("_seurat\\.rds$", "", basename(f))
  so  <- readRDS(f)
  data.frame(sample_id = sid, n_cells = ncol(so), n_genes = nrow(so),
             condition = so$condition[1], gsm_id = so$gsm_id[1],
             stringsAsFactors = FALSE)
}))

message("\n========== LOAD VERIFICATION ==========")
print(verify_df, row.names = FALSE)
message(sprintf("\nTotal cells loaded: %s",
                format(sum(verify_df$n_cells), big.mark = ",")))
message("Seurat objects saved to: ", OUT_DIR)
