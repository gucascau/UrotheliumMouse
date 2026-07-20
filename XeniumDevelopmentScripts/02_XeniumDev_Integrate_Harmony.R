################################################################################
# 02_XeniumDev_Integrate_Harmony.R
#
# Integrate the 4 Xenium developmental-kidney samples (GSE286051: W12/W92,
# male/female -- see 01_XeniumDev_BuildMetadata.R) with Harmony, then transfer
# cell type labels from the same two references used for the 12-sample UUO
# Xenium pipeline (SpatialScripts/02_Xenium_integrate_harmony.R):
#   - MKA_seurat.rds       (MouseKidneyATLAS, author_cell_type)
#   - LakesnRNA_seurat.rds (Lake IRI snRNA-seq, SubclassLevel1 + SubclassLevel3)
#
# Differs from that 12-sample script in one structural way: these 4 samples
# are flat GEO supplementary files, not a standard Xenium outs/ bundle, and
# have no cell_boundaries.parquet (no polygon segmentation) -- only
# cells.csv.gz centroids. Each sample is loaded manually (Read10X for counts,
# CreateCentroids + CreateFOV(type = "centroids") for coordinates) rather
# than via LoadXenium(). Downstream steps (normalize/PCA/Harmony/cluster/
# label transfer) are otherwise identical to that script.
#
# Input:  UsedSpatialData/Xenium/GSE286051/GSM871700[6-9]_*.{mtx,barcodes,
#         features,cells.csv}.gz
#         XeniumDevelopmentScripts/output/XeniumDev_sample_metadata.csv
#         reference/MKA_seurat.rds, reference/LakesnRNA_seurat.rds
# Output: XeniumDev_harmony_integrated_prelabeltransfer.rds
#         XeniumDev_harmony_integrated.rds
#         XeniumDev_harmony_clusters.pdf
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

options(future.globals.maxSize = 32 * 1024^3)

log_mem <- function(label) {
  g  <- gc(verbose = FALSE)
  mb <- sum(g[, 2])
  message(sprintf("  [mem] %s: %.1f GB in use", label, mb / 1024))
}

SCRIPT_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/XeniumDevelopmentScripts"
OUT_DIR    <- file.path(SCRIPT_DIR, "output")
STAGING_DIR <- file.path(OUT_DIR, "staging")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(STAGING_DIR, showWarnings = FALSE, recursive = TRUE)

REF_MKA  <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/reference/MKA_seurat.rds"
REF_LAKE <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/reference/LakesnRNA_seurat.rds"

OUT_PRELABEL <- file.path(OUT_DIR, "XeniumDev_harmony_integrated_prelabeltransfer.rds")
OUT_RDS      <- file.path(OUT_DIR, "XeniumDev_harmony_integrated.rds")
OUT_PDF      <- file.path(OUT_DIR, "XeniumDev_harmony_clusters.pdf")

# ── Load metadata ────────────────────────────────────────────────────────────
meta_csv <- file.path(OUT_DIR, "XeniumDev_sample_metadata.csv")
if (!file.exists(meta_csv)) {
  stop("Missing ", meta_csv, " -- run 01_XeniumDev_BuildMetadata.R first.")
}
sample_meta <- read.csv(meta_csv, stringsAsFactors = FALSE)
message(sprintf("==> Loaded metadata for %d samples", nrow(sample_meta)))

# ── Load one sample: counts (Read10X) + centroids-only FOV ─────────────────
load_xenium_flat_sample <- function(row) {
  message(sprintf("  Loading %s ...", row$sample_id))
  stage_dir <- file.path(STAGING_DIR, row$sample_id)
  dir.create(stage_dir, showWarnings = FALSE, recursive = TRUE)

  p <- row$file_prefix
  d <- row$matrix_dir
  file.copy(file.path(d, paste0(p, ".matrix.mtx.gz")),   file.path(stage_dir, "matrix.mtx.gz"), overwrite = TRUE)
  file.copy(file.path(d, paste0(p, ".barcodes.tsv.gz")), file.path(stage_dir, "barcodes.tsv.gz"), overwrite = TRUE)
  file.copy(file.path(d, paste0(p, ".features.tsv.gz")), file.path(stage_dir, "features.tsv.gz"), overwrite = TRUE)

  counts <- Read10X(data.dir = stage_dir, gene.column = 2)
  if (is.list(counts)) counts <- counts[["Gene Expression"]]

  cells_df <- read.csv(gzfile(file.path(d, paste0(p, ".cells.csv.gz"))), stringsAsFactors = FALSE)
  stopifnot(all(colnames(counts) %in% cells_df$cell_id))
  cells_df <- cells_df[match(colnames(counts), cells_df$cell_id), ]

  obj <- CreateSeuratObject(counts = counts, assay = "Xenium", project = row$sample_id)
  obj <- subset(obj, subset = nCount_Xenium > 10 & nFeature_Xenium > 3)

  centroids_df <- data.frame(
    x    = cells_df$x_centroid,
    y    = cells_df$y_centroid,
    cell = cells_df$cell_id
  )
  centroids_df <- centroids_df[centroids_df$cell %in% Cells(obj), ]
  cents <- CreateCentroids(centroids_df)
  fov <- CreateFOV(coords = list(centroids = cents), type = "centroids", assay = "Xenium")
  obj[["fov"]] <- fov

  obj$sample_id  <- row$sample_id
  obj$Age        <- row$Age
  obj$sex        <- row$sex
  obj$GSM        <- row$GSM

  log_mem(paste("after loading", row$sample_id))
  obj
}

# ── Step 1: Load all 4 samples ──────────────────────────────────────────────
message("==> Loading Xenium samples ...")
sample_list <- vector("list", nrow(sample_meta))
for (i in seq_len(nrow(sample_meta))) {
  sample_list[[i]] <- load_xenium_flat_sample(sample_meta[i, ])
}
names(sample_list) <- sample_meta$sample_id

cell_counts <- sapply(sample_list, ncol)
message("  Cell counts per sample:")
print(cell_counts)
message(sprintf("  Total cells: %d", sum(cell_counts)))

# ── Step 2: Merge all 4 samples ─────────────────────────────────────────────
message("==> Merging all 4 samples ...")
object <- merge(
  x            = sample_list[[1]],
  y            = sample_list[-1],
  add.cell.ids = names(sample_list),
  project      = "Xenium_Development"
)
rm(sample_list); gc()
log_mem("after merge")

DefaultAssay(object) <- "Xenium"

# ── Step 3: Normalize (matching the 12-sample UUO pipeline exactly) ────────
message("==> Normalizing ...")
object <- NormalizeData(object,
  normalization.method = "LogNormalize",
  scale.factor         = 100
)
object <- JoinLayers(object)
log_mem("after normalization")

# ── Step 4: PCA on all cells ─────────────────────────────────────────────────
message("==> Scaling and running PCA ...")
xenium_genes <- rownames(object)
message(sprintf("  %d Xenium genes used for PCA", length(xenium_genes)))

object <- ScaleData(object, features = xenium_genes)
object <- RunPCA(object,
  features       = xenium_genes,
  reduction.name = "pca",
  npcs           = 30,
  verbose        = FALSE
)
log_mem("after PCA")

# ── Step 5: Harmony batch correction ────────────────────────────────────────
message("==> Running Harmony (batch = sample_id) ...")
object <- RunHarmony(
  object         = object,
  group.by.vars  = "sample_id",
  reduction      = "pca",
  reduction.save = "harmony",
  theta          = 1,
  max_iter       = 20,
  verbose        = TRUE
)
log_mem("after Harmony")

# ── Step 6: Cluster and UMAP ─────────────────────────────────────────────────
message("==> Clustering and UMAP ...")
object <- FindNeighbors(object, reduction = "harmony", dims = 1:30)
object <- FindClusters(object, resolution = 0.3, cluster.name = "seurat_clusters")
object <- RunUMAP(object,
  reduction      = "harmony",
  dims           = 1:30,
  reduction.name = "umap"
)
log_mem("after clustering")
saveRDS(object, OUT_PRELABEL)
message("  Saved: ", OUT_PRELABEL)

# ── Reference prep helper (identical to 02_Xenium_integrate_harmony.R) ─────
prepare_reference <- function(ref_path, label) {
  message(sprintf("==> Loading %s reference ...", label))
  ref <- readRDS(ref_path)
  DefaultAssay(ref) <- "originalexp"

  feat_name <- as.character(ref[["originalexp"]]@meta.features$feature_name)
  symbol    <- sub("_ENSMUSG.*$", "", feat_name)
  keep      <- !duplicated(symbol) & symbol != ""
  ref       <- ref[keep, ]
  symbol    <- symbol[keep]

  a <- ref[["originalexp"]]
  if (nrow(a@data) > 0)       rownames(a@data)       <- symbol
  if (nrow(a@counts) > 0)     rownames(a@counts)      <- symbol
  if (nrow(a@scale.data) > 0) rownames(a@scale.data)  <- symbol
  rownames(a@meta.features) <- symbol
  ref[["originalexp"]] <- a

  shared <- intersect(xenium_genes, rownames(ref))
  message(sprintf("  Shared genes with Xenium panel: %d / %d", length(shared), length(xenium_genes)))

  ref <- ref[shared, ]
  ref <- ScaleData(ref, features = shared, verbose = FALSE)
  ref <- RunPCA(ref,
    features       = shared,
    reduction.name = "pca",
    npcs           = 30,
    verbose        = FALSE
  )
  log_mem(paste("after", label, "reference PCA"))
  list(ref = ref, shared = shared)
}

# ── Step 7: Label transfer from MKA reference ───────────────────────────────
mka        <- prepare_reference(REF_MKA, "MKA")
ref_mka    <- mka$ref
shared_mka <- mka$shared
rm(mka)

DefaultAssay(object) <- "Xenium"
message("==> FindTransferAnchors (MKA) ...")
anchors_mka <- FindTransferAnchors(
  reference            = ref_mka,
  query                = object,
  normalization.method = "LogNormalize",
  reference.reduction  = "pca",
  features             = shared_mka,
  dims                 = 1:30,
  k.anchor             = 5
)
message("==> MapQuery (MKA -> MKA_celltype) ...")
object <- MapQuery(
  anchorset           = anchors_mka,
  query               = object,
  reference           = ref_mka,
  refdata             = list(MKA_celltype = "author_cell_type"),
  reference.reduction = "pca"
)
rm(ref_mka, anchors_mka); gc()
log_mem("after MKA label transfer")

# ── Step 8: Label transfer from Lake snRNA-seq reference ───────────────────
lake        <- prepare_reference(REF_LAKE, "Lake")
ref_lake    <- lake$ref
shared_lake <- lake$shared
rm(lake)

DefaultAssay(object) <- "Xenium"
message("==> FindTransferAnchors (Lake) ...")
anchors_lake <- FindTransferAnchors(
  reference            = ref_lake,
  query                = object,
  normalization.method = "LogNormalize",
  reference.reduction  = "pca",
  features             = shared_lake,
  dims                 = 1:30,
  k.anchor             = 5
)
message("==> MapQuery (Lake -> Lake_SubclassLevel1 + Lake_SubclassLevel3) ...")
object <- MapQuery(
  anchorset           = anchors_lake,
  query               = object,
  reference           = ref_lake,
  refdata             = list(
    Lake_SubclassLevel1 = "SubclassLevel1",
    Lake_SubclassLevel3 = "SubclassLevel3"
  ),
  reference.reduction = "pca"
)
rm(ref_lake, anchors_lake); gc()
log_mem("after Lake label transfer")

# ── Step 9: Visualize and save ──────────────────────────────────────────────
message("==> Saving plots ...")
pdf(OUT_PDF, width = 16, height = 7)

p_cluster <- DimPlot(object, reduction = "umap", group.by = "seurat_clusters",
  label = TRUE, repel = TRUE, pt.size = 0.1) + ggtitle("Harmony clusters")
p_mka <- DimPlot(object, reduction = "umap", group.by = "predicted.MKA_celltype",
  label = TRUE, repel = TRUE, pt.size = 0.1) + ggtitle("MKA: author_cell_type")
p_lake1 <- DimPlot(object, reduction = "umap", group.by = "predicted.Lake_SubclassLevel1",
  label = TRUE, repel = TRUE, pt.size = 0.1) + ggtitle("Lake: SubclassLevel1 (broad)")
p_lake3 <- DimPlot(object, reduction = "umap", group.by = "predicted.Lake_SubclassLevel3",
  label = TRUE, repel = TRUE, pt.size = 0.1) + ggtitle("Lake: SubclassLevel3 (detailed)")
print(p_cluster | p_mka)
print(p_lake1 | p_lake3)

p_sample <- DimPlot(object, reduction = "umap", group.by = "sample_id",
  label = FALSE, pt.size = 0.1) + ggtitle("Sample identity")
p_age <- DimPlot(object, reduction = "umap", group.by = "Age",
  label = FALSE, pt.size = 0.1) + ggtitle("Age (W12 vs W92)")
p_sex <- DimPlot(object, reduction = "umap", group.by = "sex",
  label = FALSE, pt.size = 0.1) + ggtitle("Sex")
print(p_sample | p_age | p_sex)

for (i in seq_len(nrow(sample_meta))) {
  samp    <- sample_meta$sample_id[i]
  cells_s <- WhichCells(object, expression = sample_id == samp)
  sub_obj <- subset(object, cells = cells_s)
  tryCatch({
    print(
      ImageDimPlot(sub_obj, fov = Images(sub_obj)[1], group.by = "predicted.MKA_celltype",
                   cols = "polychrome", axes = TRUE, size = 0.3) +
        ggtitle(paste("Cell types:", samp))
    )
  }, error = function(e) {
    message("  ImageDimPlot failed for ", samp, ": ", e$message)
  })
}

dev.off()
message("  Saved: ", OUT_PDF)

message("==> Saving RDS ...")
saveRDS(object, OUT_RDS)
message("  Saved: ", OUT_RDS)

unlink(STAGING_DIR, recursive = TRUE)

message("==> Done.")
log_mem("final")
