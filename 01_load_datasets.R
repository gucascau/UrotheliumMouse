################################################################################
# Script 01: Load All Single-Cell Datasets and Create Seurat Objects
#
# Datasets overview:
# ┌────────────────────────┬───────────────┬───────────────┬────────────────────┐
# │ sample_id              │ Accession     │ Technology    │ Condition          │
# ├────────────────────────┼───────────────┼───────────────┼────────────────────┤
# │ KidneyHealthy1         │ GSE119531     │ Drop-seq      │ Healthy            │
# │ KidneyUUO1             │ GSE119531     │ snRNA-seq     │ UUO                │
# │ KidneyE18_5_1          │ GSE108291     │ 10X v2        │ Normal (E18.5)     │
# │ KidneyHealthy2         │ GSM4151577    │ 10X v2        │ Sham               │
# │ KidneyUUO2             │ GSM4151578    │ 10X v2        │ UUO day 2          │
# │ KidneyUUO3             │ GSM4151579    │ 10X v2        │ UUO day 7          │
# │ KidneyrUUO1            │ GSM4151580    │ 10X v2        │ rUUO               │
# │ KidneyHealthy3         │ GSM5333084    │ 10X v3.1      │ Sham               │
# │ KidneyUUO4             │ GSM5333085    │ 10X v3.1      │ UUO 10 days        │
# │ KidneyUUO5             │ GSM5333086    │ 10X v3.1      │ UUO 10 days        │
# │ BladderHomogenate1     │ GSM3723360    │ 10X v2        │ Bladder homogenate │
# │ BladderHomogenate2     │ GSM3723361    │ 10X v2        │ Bladder homogenate │
# │ KidneyHealthy4         │ GSM8417119    │ snRNA-seq     │ Sham               │
# │ KidneyHealthy5         │ GSM8417120    │ snRNA-seq     │ Sham               │
# │ KidneyUUO6             │ GSM8417121    │ snRNA-seq     │ UUO 4 days         │
# │ KidneyTET2UUO          │ GSM8417122    │ snRNA-seq     │ TET2KO UUO 4 days  │
# │ EmbryosE9_5ToE13_5     │ GSE119945     │ sci-RNA-seq3  │ Development        │
# │ KidneyUUO7             │ GSE264184     │ Dense matrix  │ UUO 10 days        │
# │ KidneyUUO8             │ GSE264184     │ Dense matrix  │ UUO 10 days        │
# │ KudoUUOUrothelium      │ KUDORDS       │ 10X           │ Urothelium organoid│
# │ MKA                    │ MKARDS        │ Mixed1        │ Atlas              │
# │ ChenSpatial            │ ChenRDS       │ Mixed2        │ Sex-specific       │
# │ LakesnRNA              │ LakesnRDS     │ Mixed3        │ Reference          │
# └────────────────────────┴───────────────┴───────────────┴────────────────────┘
################################################################################

# Load Bioconductor packages FIRST so they don't mask SeuratObject::Assays
# when loaded after Seurat (SingleCellExperiment -> SummarizedExperiment masks it)
library(SingleCellExperiment)
library(DropletUtils)  # for emptyDrops (unfiltered datasets)
library(Matrix)
library(data.table)
library(dplyr)
library(Seurat)        # v5+ — loaded last so SeuratObject::Assays is not masked

# Set working directory to dataset folder
DATA_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
setwd(DATA_DIR)

# Output directory for saved Seurat objects
OUT_DIR <- file.path(DATA_DIR, "seurat_objects")
dir.create(OUT_DIR, showWarnings = FALSE)

################################################################################
# SECTION 0: Master Metadata Table
# Fill in unknown fields (marked with "TODO") based on your knowledge of the data
################################################################################

sample_metadata <- data.frame(
  sample_id    = c("KidneyHealthy1", "KidneyUUO1",
                   "KidneyE18_5_1", "HumanOrg1", "HumanOrg2",
                   "KidneyHealthy2", "KidneyUUO2", "KidneyUUO3", "KidneyrUUO1",
                   "KidneyHealthy3", "KidneyUUO4", "KidneyUUO5",
                   "BladderHomogenate1", "BladderHomogenate2",
                   "KidneyHealthy4", "KidneyHealthy5", "KidneyUUO6", "KidneyTET2UUO",
                   "EmbryosE9_5ToE13_5",
                   "KidneyUUO7", "KidneyUUO8",
                   "KudoUUOUrothelium", "MKA", "ChenSpatial", "LakesnRNA"),
  gsm_id       = c("GSE119531", "GSE119531",
                   "GSE108291", "GSE108291", "GSE108291",
                   "GSM4151577", "GSM4151578", "GSM4151579", "GSM4151580",
                   "GSM5333084", "GSM5333085", "GSM5333086",
                   "GSM3723360", "GSM3723361",
                   "GSM8417119", "GSM8417120", "GSM8417121", "GSM8417122",
                   "GSE119945",
                   "GSE264184", "GSE264184",
                   "KUDORDS", "MKARDS", "ChenRDS", "LakesnRDS"),
  technology   = c("Drop-seq", "snRNA-seq",        # KidneyHealthy1=Drop-seq; KidneyUUO1=TenxNuc (10X Chromium snRNA-seq)
                   "10X_v2", "10X_v2", "10X_v2",  # GSE108291
                   "10X_v2", "10X_v2", "10X_v2", "10X_v2",   # GSM4151577-80
                   "10X_v3.1", "10X_v3.1", "10X_v3.1",       # GSM5333084-86: Chromium Next GEM 3' v3.1
                   "10X_v2", "10X_v2",             # GSM3723360-61
                   "snRNA-seq", "snRNA-seq", "snRNA-seq", "snRNA-seq",  # GSM8417119-22
                   "sci-RNA-seq3",                 # GSE119945: confirmed sci-RNA-seq3
                   "scRNA-seq", "scRNA-seq",        # GSE264184: 10X Chromium scRNA-seq (TODO: verify sc vs sn from paper PMID_35941120)
                   "10X", "TODO", "TODO", "snRNA-seq"),
  species      = c("mouse", "mouse",
                   "mouse", "human", "human",
                   "mouse", "mouse", "mouse", "mouse",
                   "mouse", "mouse", "mouse",
                   "mouse", "mouse",
                   "mouse", "mouse", "mouse", "mouse",
                   "mouse",
                   "mouse", "mouse",
                   "mouse", "mouse", "mouse", "mouse"),
  tissue       = c("kidney", "kidney",
                   "kidney_E18.5", "kidney_organoid", "kidney_organoid",
                   "kidney", "kidney", "kidney", "kidney",
                   "kidney", "kidney", "kidney",
                   "bladder", "bladder",
                   "kidney", "kidney", "kidney", "kidney",
                   "AllOrganoids",
                   "kidney", "kidney",
                   "kidney", "kidney", "kidney", "kidney"),
  condition    = c("Healthy", "UUO",
                   "Normal", "Normal", "Normal",
                   "Sham", "UUO_day2", "UUO_day7", "rUUO",
                   "Sham", "UUO_10days", "UUO_10days",
                   "Bladder_Homogenate", "Bladder_Homogenate",
                   "Sham", "Sham", "UUO_4days", "TET2KO_UUO_4days",
                   "Development",
                   "UUO_10days", "UUO_10days",
                   "Urothelium_Organoid", "Atlas_1", "Atlas_2", "Atlas_3"),
  paper        = c("PMID_30510133", "PMID_30510133",
                   "PMID_31118232", "PMID_31118232", "PMID_31118232",
                   "PMID_32978267;PMID_36509292", "PMID_32978267;PMID_36509292", "PMID_32978267;PMID_36509292", "PMID_32978267;PMID_36509292",
                   "PMID_39414946", "PMID_39414946", "PMID_39414946",
                   "PMID_31462402", "PMID_31462402",
                   "PMID_39511169", "PMID_39511169", "PMID_39511169", "PMID_39511169",
                   "PMID_30787437",
                   "PMID_35941120", "PMID_35941120",
                   "Kudo_et_al", "PMID_37275529;PMID_41071604", "PMID_40259083", "PMID_41040183"),
  has_annotation = c(TRUE, TRUE,
                     FALSE, FALSE, FALSE,
                     FALSE, FALSE, FALSE, FALSE,
                     FALSE, FALSE, FALSE,
                     FALSE, FALSE,
                     FALSE, FALSE, FALSE, FALSE,
                     TRUE,
                     FALSE, FALSE,
                     TRUE, TRUE, TRUE, TRUE),
  stringsAsFactors = FALSE
)


################################################################################
# SECTION 1: Helper Functions
################################################################################

# Helper: Load standard 10X MTX (barcodes/genes/matrix)
load_10x_mtx <- function(barcode_file, gene_file, matrix_file,
                          gene_col = 2,  # column in gene file to use as rownames
                          min_cells = 3, min_features = 200) {
  bc_con <- file(barcode_file, "rb"); on.exit(close(bc_con), add = TRUE)
  barcodes <- readLines(gzcon(bc_con))

  genes_df <- read.table(gzfile(gene_file), header = FALSE, sep = "\t",
                          stringsAsFactors = FALSE)
  gene_names <- make.unique(genes_df[, gene_col])

  mat_con <- file(matrix_file, "rb"); on.exit(close(mat_con), add = TRUE)
  mat <- readMM(gzcon(mat_con))
  rownames(mat) <- gene_names
  colnames(mat) <- barcodes
  return(mat)
}

# Helper: Filter unfiltered 10X matrix using emptyDrops (DropletUtils)
filter_unfiltered_10x <- function(mat, fdr_threshold = 0.01, lower = 100) {
  message("Running emptyDrops for ambient RNA filtering...")
  set.seed(42)
  e_out <- emptyDrops(mat, lower = lower)
  is_cell <- !is.na(e_out$FDR) & e_out$FDR <= fdr_threshold
  message(sprintf("  Retained %d cells (FDR <= %g)", sum(is_cell), fdr_threshold))
  return(mat[, is_cell])
}

# Helper: Create Seurat object and add metadata
make_seurat <- function(mat, sample_id, meta_row, min_cells = 3, min_features = 200) {
  # Seurat v5 requires rownames on the count matrix; fail early with a clear message
  if (is.null(rownames(mat)) || any(rownames(mat) == ""))
    stop(sprintf("[make_seurat] Matrix for '%s' has missing rownames — check loading code.", sample_id))
  so <- CreateSeuratObject(
    counts      = mat,
    project     = sample_id,
    min.cells   = min_cells,
    min.features = min_features
  )
  so$sample_id    <- sample_id
  so$gsm_id       <- meta_row$gsm_id
  so$technology   <- meta_row$technology
  so$species      <- meta_row$species
  so$tissue       <- meta_row$tissue
  so$condition    <- meta_row$condition
  so$paper        <- meta_row$paper
  so$has_annotation <- meta_row$has_annotation
  return(so)
}

get_meta <- function(sid) sample_metadata[sample_metadata$sample_id == sid, ]

# Helper: Create Seurat object for pre-normalized data (e.g. GSE264184)
# Puts values in the data layer, not counts, so NormalizeData must be skipped.
# Sets $skip_normalization = TRUE as a flag for script 02.
make_prenorm_seurat <- function(mat, sample_id, meta_row, min_cells = 3, min_features = 200) {
  cells_ok <- Matrix::colSums(mat > 0) >= min_features
  genes_ok  <- Matrix::rowSums(mat[, cells_ok] > 0) >= min_cells
  mat_filt  <- mat[genes_ok, cells_ok]
  assay_obj <- CreateAssay5Object(data = mat_filt)
  so <- CreateSeuratObject(counts = assay_obj, project = sample_id)
  so$sample_id          <- sample_id
  so$gsm_id             <- meta_row$gsm_id
  so$technology         <- meta_row$technology
  so$species            <- meta_row$species
  so$tissue             <- meta_row$tissue
  so$condition          <- meta_row$condition
  so$paper              <- meta_row$paper
  so$has_annotation     <- meta_row$has_annotation
  so$skip_normalization <- TRUE
  return(so)
}


################################################################################
# SECTION 2: GSE119531 – Drop-seq (mouse kidney Healthy + UUO)
# Ref: Wu H, Kirita Y, Donnelly EL, Humphreys BD. Advantages of Single-Nucleus over Single-Cell RNA Sequencing of Adult Kidney: Rare Cell Types and Novel Cell States Revealed in Fibrosis. J Am Soc Nephrol 2019 Jan;30(1):23-32. PMID: 30510133
# Format: dense tab-delimited matrix (genes x cells), separate annotation file
# Known cell types: PT(S1-S2), PT(S3), DL/tAL, LOH, DCT, CD, Endo, Fib, Immune
################################################################################

load_gse119531 <- function(condition = "Healthy") {
  if (condition == "Healthy") {
    dge_file  <- "GSE119531_Healthy.combined.dge.txt"
    ann_file  <- "GSE119531_Healthy.combined.cell.annotation.txt"
    sample_id <- "KidneyHealthy1"
  } else {
    dge_file  <- "GSE119531_UUO.dge.txt"
    ann_file  <- "GSE119531_UUO.cell.annotation.txt"
    sample_id <- "KidneyUUO1"
  }

  message(sprintf("Loading GSE119531 %s...", condition))
  # read.table with row.names=1 correctly handles the Drop-seq DGE format where
  # the header row has an empty first field (tab-leading), which fread misparses.
  # check.names=FALSE preserves barcode strings exactly as-is.
  mat_df <- read.table(dge_file, header = TRUE, sep = "\t", row.names = 1,
                       stringsAsFactors = FALSE, check.names = FALSE)
  genes <- make.unique(rownames(mat_df))
  rownames(mat_df) <- genes

  message(sprintf("  Dims: %d genes x %d cells", nrow(mat_df), ncol(mat_df)))
  message(sprintf("  First 3 genes   : %s", paste(head(genes, 3), collapse = ", ")))
  message(sprintf("  First 3 barcodes: %s", paste(head(colnames(mat_df), 3), collapse = ", ")))
  message(sprintf("  Empty gene names: %s | NA gene names: %s",
                  any(genes == ""), any(is.na(genes))))

  mat_sparse <- as(as.matrix(mat_df), "dgCMatrix")

  ann <- read.table(ann_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  # Diagnostics: print column names of annotation file and first few barcodes
  message(sprintf("  Ann columns: %s", paste(colnames(ann), collapse = ", ")))
  message(sprintf("  First 3 ann barcodes : %s", paste(head(ann[[1]], 3), collapse = ", ")))
  message(sprintf("  First 3 mat barcodes : %s", paste(head(colnames(mat_sparse), 3), collapse = ", ")))
  message(sprintf("  Matrix dims: %d genes x %d cells", nrow(mat_sparse), ncol(mat_sparse)))

  # Use first column of ann as barcode column (handles different column names)
  barcode_col  <- ann[[1]]

  # The matrix prefix (e.g. "UUO1_") may differ from the annotation prefix
  # (e.g. "UUO_"). Strip everything up to and including the first "_" from both
  # sides, then match on the bare 16-bp barcode sequence.
  mat_suffix <- sub("^[^_]+_", "", colnames(mat_sparse))
  ann_suffix <- sub("^[^_]+_", "", barcode_col)
  common_suffix <- intersect(mat_suffix, ann_suffix)
  message(sprintf("  Common cells: %d", length(common_suffix)))

  if (length(common_suffix) == 0)
    stop(sprintf("[load_gse119531] No matching barcodes for %s after prefix stripping.", sample_id))

  mat_idx <- match(common_suffix, mat_suffix)
  ann_idx <- match(common_suffix, ann_suffix)

  gene_names   <- rownames(mat_sparse)
  mat_sparse   <- mat_sparse[, mat_idx, drop = FALSE]
  rownames(mat_sparse) <- gene_names
  # Rename columns to match annotation barcodes exactly
  colnames(mat_sparse) <- barcode_col[ann_idx]
  ann <- ann[ann_idx, ]

  so <- make_seurat(mat_sparse, sample_id, get_meta(sample_id))
  so$cell_type_original <- ann$CellType[match(colnames(so), barcode_col)]
  return(so)
}

so_healthy <- load_gse119531("Healthy")
so_uuo     <- load_gse119531("UUO")


################################################################################
# SECTION 3: GSE108291 – 10X v2 (E18.5 kidneys + human organoids)
# Ref: Combes AN, Phipson B, Lawlor KT, Dorison A et al. Single cell analysis of the developing mouse kidney provides deeper insight into marker gene expression and ligand-receptor crosstalk. Development 2019 Jun 12;146(12). PMID: 31118232
# kid: Mouse kidney (ENSMUSG IDs) – integrate with mouse datasets
# org / org4: Human kidney organoids (ENSG IDs) – requires ortholog mapping
#             Set INCLUDE_HUMAN_ORGANOIDS = TRUE to include
################################################################################

INCLUDE_HUMAN_ORGANOIDS <- FALSE  # <-- set TRUE if you want cross-species integration

load_gse108291_kid <- function() {
  message("Loading GSE108291 kidney (mouse)...")
  mat <- load_10x_mtx(
    barcode_file = "GSE108291_kid_barcodes.tsv.gz",
    gene_file    = "GSE108291_kid_genes.tsv.gz",
    matrix_file  = "GSE108291_kid_matrix.mtx.gz",
    gene_col     = 2
  )
  # NOTE: This is an unfiltered matrix (2.2M barcodes). Filter with emptyDrops.
  mat_filt <- filter_unfiltered_10x(mat)
  make_seurat(mat_filt, "KidneyE18_5_1", get_meta("KidneyE18_5_1"))
}

so_kid <- load_gse108291_kid()

if (INCLUDE_HUMAN_ORGANOIDS) {
  # Human organoids – map human genes to mouse orthologs before integration
  # Requires biomaRt:
  #   library(biomaRt)
  #   human  <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  #   mouse  <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  #   orthologs <- getLDS(attributes  = "hgnc_symbol", mart = human,
  #                       attributesL = "mgi_symbol", martL = mouse)
  message("Loading GSE108291 organoids (human – ortholog mapping required)...")
  mat_org <- load_10x_mtx("GSE108291_org_barcodes.tsv.gz",
                           "GSE108291_org_genes.tsv.gz",
                           "GSE108291_org_matrix.mtx.gz", gene_col = 2)
  mat_org4 <- load_10x_mtx("GSE108291_org4_barcodes.tsv.gz",
                            "GSE108291_org4_genes.tsv.gz",
                            "GSE108291_org4_matrix.mtx.gz", gene_col = 2)
  # After ortholog mapping, create Seurat objects:
  # so_org  <- make_seurat(mat_org_mapped,  "HumanOrg1", get_meta("HumanOrg1"))
  # so_org4 <- make_seurat(mat_org4_mapped, "HumanOrg2", get_meta("HumanOrg2"))
}


################################################################################
# SECTION 4: GSM4151577–4151580 – 10X v2 (Kirita et al. 2020, mouse kidney UUO)
# Ref: Conway BR, O'Sullivan ED, Cairns C, O'Sullivan J et al. Kidney Single-Cell Atlas Reveals Myeloid Heterogeneity in Progression and Regression of Kidney Disease. J Am Soc Nephrol 2020 Dec;31(12):2833-2854. PMID: 32978267
# O'Sullivan ED, Mylonas KJ, Bell R, Carvalho C et al. Single-cell analysis of senescent epithelia reveals targetable mechanisms promoting fibrosis. JCI Insight 2022 Nov 22;7(22). PMID: 36509292
# Samples: sham, uuo2, uuo7, ruuo  (filtered, ~3–6k cells each)
################################################################################

kirita_samples <- list(
  KidneyHealthy2 = list(b = "GSM4151577_sham_barcodes.tsv.gz",
                         g = "GSM4151577_sham_genes.tsv.gz",
                         m = "GSM4151577_sham_matrix.mtx.gz"),
  KidneyUUO2     = list(b = "GSM4151578_uuo2_barcodes.tsv.gz",
                         g = "GSM4151578_uuo2_genes.tsv.gz",
                         m = "GSM4151578_uuo2_matrix.mtx.gz"),
  KidneyUUO3     = list(b = "GSM4151579_uuo7_barcodes.tsv.gz",
                         g = "GSM4151579_uuo7_genes.tsv.gz",
                         m = "GSM4151579_uuo7_matrix.mtx.gz"),
  KidneyrUUO1    = list(b = "GSM4151580_ruuo_barcodes.tsv.gz",
                         g = "GSM4151580_ruuo_genes.tsv.gz",
                         m = "GSM4151580_ruuo_matrix.mtx.gz")
)

so_kirita <- lapply(names(kirita_samples), function(sid) {
  message(sprintf("Loading Kirita %s...", sid))
  f <- kirita_samples[[sid]]
  mat <- load_10x_mtx(f$b, f$g, f$m, gene_col = 2)
  make_seurat(mat, sid, get_meta(sid))
})
names(so_kirita) <- names(kirita_samples)


################################################################################
# SECTION 5: GSM5333084–5333086 – 10X v3 (Cell Ranger 4, CS/CU samples)
# Samples: CS1KKY (ctrl sham), CU1KKY (CUU rep1), CU2KKY (CUU rep2)
# NOTE: Unfiltered (6.8M barcodes each) – use emptyDrops
################################################################################

cukky_samples <- list(
  KidneyHealthy3 = list(b = "GSM5333084_CS1KKY_barcodes.tsv.gz",
                          g = "GSM5333084_CS1KKY_features.tsv.gz",
                          m = "GSM5333084_CS1KKY_matrix.mtx.gz"),
  KidneyUUO4     = list(b = "GSM5333085_CU1KKY_barcodes.tsv.gz",
                          g = "GSM5333085_CU1KKY_features.tsv.gz",
                          m = "GSM5333085_CU1KKY_matrix.mtx.gz"),
  KidneyUUO5     = list(b = "GSM5333086_CU2KKY_barcodes.tsv.gz",
                          g = "GSM5333086_CU2KKY_features.tsv.gz",
                          m = "GSM5333086_CU2KKY_matrix.mtx.gz")
)

so_cukky <- lapply(names(cukky_samples), function(sid) {
  message(sprintf("Loading CS/CU sample %s (unfiltered – running emptyDrops)...", sid))
  f <- cukky_samples[[sid]]
  mat <- load_10x_mtx(f$b, f$g, f$m, gene_col = 2)
  mat_filt <- filter_unfiltered_10x(mat)
  make_seurat(mat_filt, sid, get_meta(sid))
})
names(so_cukky) <- names(cukky_samples)


################################################################################
# SECTION 6: GSM3723360–3723361 – 10X v2 (Sample4, Sample5)
# Ref: Yu Z, Liao J, Chen Y, Zou C et al. Single-Cell Transcriptomic Map of the Human and Mouse Bladders. J Am Soc Nephrol 2019 Nov;30(11):2159-2176. PMID: 31462402
################################################################################

so_s4 <- {
  message("Loading GSM3723360 BladderHomogenate1...")
  mat <- load_10x_mtx("GSM3723360_Sample4_barcodes.tsv.gz",
                       "GSM3723360_Sample4_genes.tsv.gz",
                       "GSM3723360_Sample4_matrix.mtx.gz", gene_col = 2)
  make_seurat(mat, "BladderHomogenate1", get_meta("BladderHomogenate1"))
}

so_s5 <- {
  message("Loading GSM3723361 BladderHomogenate2...")
  mat <- load_10x_mtx("GSM3723361_Sample5_barcodes.tsv.gz",
                       "GSM3723361_Sample5_genes.tsv.gz",
                       "GSM3723361_Sample5_matrix.mtx.gz", gene_col = 2)
  make_seurat(mat, "BladderHomogenate2", get_meta("BladderHomogenate2"))
}


################################################################################
# SECTION 7: GSM8417119–8417122 – snRNA-seq (10X Chromium, NOT Multiome)
# Ref: Liang, X., Liu, H., Hu, H. et al. TET2 germline variants promote kidney disease by impairing DNA repair and activating cytosolic nucleotide sensors. Nat Commun 15, 9621 (2024).
# GEO confirms: single-nucleus RNA-seq using 10X Chromium; Illumina NovaSeq 6000
# Features files use standard 3-column format (all rows = "Gene Expression")
# load_multiome_rna() works correctly here: rna_idx selects all features
################################################################################

load_multiome_rna <- function(barcode_f, feature_f, matrix_f, sample_id) {
  message(sprintf("Loading 10X Multiome RNA for %s...", sample_id))
  bc_con <- file(barcode_f, "rb"); on.exit(close(bc_con), add = TRUE)
  barcodes  <- readLines(gzcon(bc_con))

  feat_df   <- read.table(gzfile(feature_f), header = FALSE, sep = "\t",
                           stringsAsFactors = FALSE)
  # Feature types are in column 3
  rna_idx   <- which(feat_df[, 3] == "Gene Expression")
  gene_names <- make.unique(feat_df[rna_idx, 2])

  mat_con <- file(matrix_f, "rb"); on.exit(close(mat_con), add = TRUE)
  mat_full  <- readMM(gzcon(mat_con))
  colnames(mat_full) <- barcodes
  mat_rna   <- mat_full[rna_idx, ]
  rownames(mat_rna) <- gene_names

  make_seurat(mat_rna, sample_id, get_meta(sample_id))
}

multiome_samples <- list(
  KidneyHealthy4 = c("GSM8417119_sample3_barcodes.tsv.gz",
                      "GSM8417119_sample3_features.tsv.gz",
                      "GSM8417119_sample3_matrix.mtx.gz"),
  KidneyHealthy5 = c("GSM8417120_sample4_barcodes.tsv.gz",
                      "GSM8417120_sample4_features.tsv.gz",
                      "GSM8417120_sample4_matrix.mtx.gz"),
  KidneyUUO6     = c("GSM8417121_sample5_barcodes.tsv.gz",
                      "GSM8417121_sample5_features.tsv.gz",
                      "GSM8417121_sample5_matrix.mtx.gz"),
  KidneyTET2UUO  = c("GSM8417122_sample6_barcodes.tsv.gz",
                      "GSM8417122_sample6_features.tsv.gz",
                      "GSM8417122_sample6_matrix.mtx.gz")
)

so_multiome <- lapply(names(multiome_samples), function(sid) {
  f <- multiome_samples[[sid]]
  load_multiome_rna(f[1], f[2], f[3], sid)
})
names(so_multiome) <- names(multiome_samples)


################################################################################
# SECTION 8: GSE264184 – Dense matrix (genes x cells, D1_ / D2_ prefix)
# Ref: Wang Y, Li Y, Chen Z, Yuan Y et al. GSDMD-dependent neutrophil extracellular traps promote macrophage-to-myofibroblast transition and renal fibrosis in obstructive nephropathy. Cell Death Dis 2022 Aug 8;13(8):693. PMID: 35941120
# NOTE: Contains normalized (non-integer) values – treated as log-normalized
################################################################################

message("Loading GSE264184...")
du_mat_raw  <- fread("GSE264184_DU-matrix-data.txt.gz", header = TRUE)
du_genes    <- make.unique(as.character(du_mat_raw[[1]]))
du_df       <- as.data.frame(du_mat_raw[, -1, with = FALSE])
rownames(du_df) <- du_genes
du_sparse   <- as(as.matrix(du_df), "dgCMatrix")

# Split D1 and D2 samples by barcode prefix
d1_cells <- grep("^D1_", colnames(du_sparse), value = TRUE)
d2_cells <- grep("^D2_", colnames(du_sparse), value = TRUE)

# Values are pre-normalized floats: use make_prenorm_seurat so data goes into the
# data layer (not counts). Script 02 must skip NormalizeData when $skip_normalization == TRUE.
so_du1 <- make_prenorm_seurat(du_sparse[, d1_cells], "KidneyUUO7", get_meta("KidneyUUO7"))
so_du2 <- make_prenorm_seurat(du_sparse[, d2_cells], "KidneyUUO8", get_meta("KidneyUUO8"))


################################################################################
# SECTION 9: GSE119945 – sci-RNA-seq (Cao et al. 2019, mouse kidney development)
# Format: MatrixMarket + cell annotation CSV + gene annotation CSV
# Cao J, Spielmann M, Qiu X, Huang X et al. The single-cell transcriptional landscape of mammalian organogenesis. Nature 2019 Feb;566(7745):496-502. PMID: 30787437
# ~2M barcodes; filter to annotated (non-doublet) cells only
################################################################################

message("Loading GSE119945 (sci-RNA-seq)...")
cell_ann  <- fread("GSE119945_cell_annotate.csv")
gene_ann  <- fread("GSE119945_gene_annotate.csv")

# Filter: remove doublets and unannotated cells
cell_filt <- cell_ann[detected_doublet == FALSE & !is.na(Main_Cluster)]
message(sprintf("  Retaining %d / %d cells after doublet/annotation filter",
                nrow(cell_filt), nrow(cell_ann)))

mat_sci   <- readMM("GSE119945_gene_count.txt")
# Dimensions: genes x cells (26,183 genes x 2,058,652 barcodes)
rownames(mat_sci) <- make.unique(gene_ann$gene_short_name)
colnames(mat_sci) <- cell_ann$sample   # 'sample' column holds cell barcodes

mat_sci_filt <- mat_sci[, cell_filt$sample]

so_sci <- make_seurat(mat_sci_filt, "EmbryosE9_5ToE13_5", get_meta("EmbryosE9_5ToE13_5"))
# Use match() so annotations stay correct if CreateSeuratObject dropped any cells
# (e.g. cells that failed the min.features = 200 threshold)
so_sci$cell_type_original <- cell_filt$Sub_Cluster[match(colnames(so_sci), cell_filt$sample)]
so_sci$main_cluster       <- cell_filt$Main_Cluster[match(colnames(so_sci), cell_filt$sample)]
so_sci$sex                <- cell_filt$sex[match(colnames(so_sci), cell_filt$sample)]
so_sci$dev_day            <- cell_filt$day[match(colnames(so_sci), cell_filt$sample)]


################################################################################
# SECTION 10: Pre-processed RDS Files
# Three zellkonverter-converted SingleCellExperiment objects + one Seurat RDS
################################################################################

library(zellkonverter)

## 10a. KudoUUO – already merged and Harmony-integrated Seurat object
## Our own daaset that integrated the original Kudo UUO data (GSE37275529) with our new data (PMID_41071604).
message("Loading KudoUUO RDS (Seurat)...")
so_kudo <- readRDS("KudoUUOComparision_Scrna_Merged_RunHarmony.rds")
so_kudo$sample_id    <- "KudoUUOUrothelium"
so_kudo$gsm_id       <- get_meta("KudoUUOUrothelium")$gsm_id
so_kudo$technology   <- get_meta("KudoUUOUrothelium")$technology
so_kudo$condition    <- get_meta("KudoUUOUrothelium")$condition
so_kudo$paper        <- get_meta("KudoUUOUrothelium")$paper
so_kudo$tissue          <- "kidney"
so_kudo$species         <- "mouse"
so_kudo$has_annotation  <- get_meta("KudoUUOUrothelium")$has_annotation

## 10b. Mouse Kidney Atlas (SCE -> Seurat)
## 
message("Loading MouseKidneyAtlas RDS (SCE -> Seurat)...")
# Ref: Novella-Rausell, C., Grudniewska, M., Peters, D.J.M., Mahfouz, A. A comprehensive mouse kidney atlas enables rare cell population characterization and robust marker discovery. iScience, 26(6): 106877 (2023). Yasinoglu, S.A. and Novella-Rausell, C., et al. Spatial Transcriptomics Reveals Injured Cells, Signature Genes, and Communication Patterns in the Cyst Microenvironment of Polycystic Kidney Disease. Journal of the American Society of Nephrology, 2025.

so_mka   <- readRDS("MouseKidneyATLAS_MKA_updated_zellkonvertedConverted.rds")
so_mka$sample_id  <- "MKA"
so_mka$gsm_id     <- get_meta("MKA")$gsm_id
so_mka$technology <- get_meta("MKA")$technology
so_mka$condition  <- get_meta("MKA")$condition
so_mka$paper      <- get_meta("MKA")$paper
so_mka$tissue          <- "kidney"
so_mka$species         <- "mouse"
so_mka$has_annotation  <- get_meta("MKA")$has_annotation

## 10c. Chen 2025 NatGenet – Multi-omic Spatial (SCE -> Seurat)
# Ref: Chen S, Liu R, Mo CK, Wendl MC, Houston A, Lal P, Zhao Y, Caravan W, Shinkle AT, Abedin-Do A, Naser Al Deen N, Sato K, Li X, Targino da Costa ALN, Li Y, Karpova A, Herndon JM, Artyomov MN, Rubin JB, Jain S, Li X, Stewart SA, Ding L, Chen F. Multi-omic and spatial analysis of mouse kidneys highlights sex-specific differences in gene regulation across the lifespan. Nat Genet. 2025 May;57(5):1213-1227. doi: 10.1038/s41588-025-02161-x. Epub 2025 Apr 21. PMID: 40259083; PMCID: PMC12081296.
message("Loading Chen2025 RDS (SCE -> Seurat)...")
so_chen  <- readRDS("MultiOmicSpatialMouseKidney_SexSpecific_lifespan_Chen2025NatGenet_zellkonvertedConverted.rds")
#so_chen   <- as.Seurat(sce_chen, counts = "counts", data = "logcounts")
so_chen$sample_id  <- "ChenSpatial"
so_chen$gsm_id     <- get_meta("ChenSpatial")$gsm_id
so_chen$technology <- get_meta("ChenSpatial")$technology
so_chen$condition  <- get_meta("ChenSpatial")$condition
so_chen$paper      <- get_meta("ChenSpatial")$paper
so_chen$tissue          <- "kidney"
so_chen$species         <- "mouse"
so_chen$has_annotation  <- get_meta("ChenSpatial")$has_annotation

## 10d. Lake 2025 bioRxiv – snRNA-seq (SCE -> Seurat)
## Ref: Lake BB, Melo Ferreira R, Hansen J, Menon R, Basta J, Thiessen Philbrook H, Reinert S, Fallegger R, Lagwankar AK, Chen X, Maity S, Djambazova KV, Gorman BL, Lucarelli N, Gisch DL, Schmidt IM, Nair V, Alakwaa F, Kefaloyianni E, Zhang B, Knoten AL, Kaushal M, Otto EA, Farrow MA, Diep D, Velickovic D, Sabo AR, Cole E, Tamayo I, Tanevski J, Conklin KY, Sealfon RSG, He Y, Brennan M, Robbins L, Cheng YH, Bitzer M, Surapaneni A, Menez S, Kharchenko PV, Alpers CE, Balis UGJ, Barisoni L, de Boer IH, Demeke D, Fogo AB, Henderson JM, Herlitz L, Moeckel GW, Randhawa PS, Rosenberg AZ, Schaub JA, Setty S, Brosius FC, Caramori ML, Coca SG, Figenshau RS, Kim EH, Kiryluk K, Lash JP, Miller RT, O'Toole JF, Palevsky PM, Rhee EP, Ricardo AC, Rosas SE, Roy-Chaudhury P, Sarwal MM, Sedor JR, Toto RD, Turkmen A, Waikar SS, Williams JC, Wilson FP, Woodle ES, Macosko EZ, Saez-Rodriguez J, Dagher PC, Grams ME, Bjornstad P, El-Achkar TM, Troyanskaya OG, Bonevich N, Sarder P, Kumar S, Anderton CR, Spraggins JM, Sharma K, Rauchman M, Himmelfarb J, Gaut JP, Precision Medicine Project K, Zhang K, Iyengar R, Kretzler M, Hodgin JB, Parikh CR, Eadon MT, Jain S. Cellular and Spatial Drivers of Unresolved Injury and Functional Decline in the Human Kidney. bioRxiv [Preprint]. 2025 Nov 24:2025.09.26.678707. doi: 10.1101/2025.09.26.678707. PMID: 41040183; PMCID: PMC12485684.
message("Loading Lake2025 RDS (SCE -> Seurat)...")
so_lake  <- readRDS("mouse_kidney_snRNAseq_Lake2025_bioRxiv_V2_zellkonvertedConverted.rds")
# so_lake   <- as.Seurat(sce_lake, counts = "counts", data = "logcounts")
so_lake$sample_id  <- "LakesnRNA"
so_lake$gsm_id     <- get_meta("LakesnRNA")$gsm_id
so_lake$technology <- get_meta("LakesnRNA")$technology
so_lake$condition  <- get_meta("LakesnRNA")$condition
so_lake$paper      <- get_meta("LakesnRNA")$paper
so_lake$tissue          <- "kidney"
so_lake$species         <- "mouse"
so_lake$has_annotation  <- get_meta("LakesnRNA")$has_annotation


################################################################################
# SECTION 11: Collect all mouse Seurat objects into a named list
################################################################################

all_seurat_objects <- c(
  # Drop-seq
  list(KidneyHealthy1 = so_healthy, KidneyUUO1 = so_uuo),
  # 10X kidney E18.5 (mouse only)
  list(KidneyE18_5_1 = so_kid),
  # Kirita UUO (KidneyHealthy2/KidneyUUO2/KidneyUUO3/KidneyrUUO1)
  so_kirita,
  # CS/CU (KidneyHealthy3/KidneyUUO4/KidneyUUO5)
  so_cukky,
  # GSM3723360-61 bladder
  list(BladderHomogenate1 = so_s4, BladderHomogenate2 = so_s5),
  # Multiome RNA (KidneyHealthy4/5, KidneyUUO6, KidneyTET2UUO)
  so_multiome,
  # GSE264184 pre-normalized
  list(KidneyUUO7 = so_du1, KidneyUUO8 = so_du2),
  # sci-RNA-seq development
  list(EmbryosE9_5ToE13_5 = so_sci),
  # Pre-processed atlases/references
  list(KudoUUOUrothelium = so_kudo, MKA = so_mka, ChenSpatial = so_chen, LakesnRNA = so_lake)
)

message(sprintf("Total datasets loaded: %d", length(all_seurat_objects)))

# Save individual Seurat objects
for (nm in names(all_seurat_objects)) {
  saveRDS(all_seurat_objects[[nm]],
          file = file.path(OUT_DIR, paste0(nm, "_seurat.rds")))
}

################################################################################
# VERIFICATION: print and save a metadata summary table to confirm all fields
# are correctly assigned. Check this in the log after the job completes.
################################################################################

meta_cols <- c("sample_id", "gsm_id", "technology", "species",
               "tissue", "condition", "paper", "has_annotation")

verify_df <- do.call(rbind, lapply(names(all_seurat_objects), function(nm) {
  so <- all_seurat_objects[[nm]]
  data.frame(
    sample_id      = nm,
    gsm_id         = so$gsm_id[1],
    technology     = so$technology[1],
    species        = so$species[1],
    tissue         = so$tissue[1],
    condition      = so$condition[1],
    paper          = so$paper[1],
    has_annotation = so$has_annotation[1],
    n_cells        = ncol(so),
    n_genes        = nrow(so),
    stringsAsFactors = FALSE
  )
}))

# Flag any remaining TODOs
todo_fields <- which(verify_df == "TODO", arr.ind = TRUE)
if (nrow(todo_fields) > 0) {
  warning("The following metadata fields are still set to 'TODO':")
  for (i in seq_len(nrow(todo_fields))) {
    warning(sprintf("  sample: %-15s | column: %s",
                    verify_df$sample_id[todo_fields[i, 1]],
                    colnames(verify_df)[todo_fields[i, 2]]))
  }
}

# Print to log
message("\n========== METADATA VERIFICATION ==========")
print(verify_df, row.names = FALSE)
message(sprintf("\nTotal cells across all datasets: %s",
                format(sum(verify_df$n_cells), big.mark = ",")))

# Save as CSV for easy inspection
write.csv(verify_df,
          file = file.path(OUT_DIR, "metadata_verification.csv"),
          row.names = FALSE)
message(sprintf("Verification table saved to: %s/metadata_verification.csv", OUT_DIR))
message("Done. Individual Seurat objects saved to: ", OUT_DIR)
