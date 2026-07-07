################################################################################
# 02_Xenium_integrate_harmony.R
#
# Integrate 12 Xenium mouse kidney UUO time-course samples with Harmony,
# then transfer cell type labels from two references separately:
#   - MKA_seurat.rds       (MouseKidneyATLAS, author_cell_type)
#   - LakesnRNA_seurat.rds (Lake IRI snRNA-seq, SubclassLevel1 + SubclassLevel3)
#
# Samples:
#   GSM8325615  ShamL    GSM8325616  ShamR
#   GSM8325617  Hour4L   GSM8325618  Hour4R
#   GSM8325619  Hour12L  GSM8325620  Hour12R
#   GSM8325621  Day2L    GSM8325622  Day2R
#   GSM8325623  Day14L   GSM8325624  Day14R
#   GSM8325625  Week6L   GSM8325626  Week6R
#
# Pipeline (following published IRI Xenium workflow):
#   1.  LoadXenium() per sample + QC filter (nCount > 10, nFeature > 3)
#   2.  Merge all 12 samples
#   3.  NormalizeData (scale.factor = 100) — simple log-norm, no SCTransform
#   4.  ScaleData + RunPCA on all ~300 Xenium genes (no sketching)
#   5.  RunHarmony (batch = sample_id) on all cells
#   6.  FindNeighbors + FindClusters + RunUMAP
#   7.  Label transfer from MKA reference  -> MKA_celltype
#   8.  Label transfer from Lake reference -> Lake_SubclassLevel1, Lake_SubclassLevel3
#   9.  Save integrated object + plots
################################################################################

library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)
library(dplyr)

options(future.globals.maxSize = 32 * 1024^3)

log_mem <- function(label) {
  g  <- gc(verbose = FALSE)
  mb <- sum(g[, 2])
  message(sprintf("  [mem] %s: %.1f GB in use", label, mb / 1024))
}

# ── Paths ──────────────────────────────────────────────────────────────────────
XEN_BASE <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSpatialData/Xenium"
REF_MKA  <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/reference/MKA_seurat.rds"
REF_LAKE <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/reference/LakesnRNA_seurat.rds"

OUT_DIR  <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/SpatialScripts/output"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OUT_RDS  <- file.path(OUT_DIR, "Xenium_harmony_integrated.rds")
OUT_PDF  <- file.path(OUT_DIR, "Xenium_harmony_clusters.pdf")

# ── Sample manifest ────────────────────────────────────────────────────────────
sample_meta <- data.frame(
  gsm       = c("GSM8325615","GSM8325616","GSM8325617","GSM8325618",
                "GSM8325619","GSM8325620","GSM8325621","GSM8325622",
                "GSM8325623","GSM8325624","GSM8325625","GSM8325626"),
  sample_id = c("ShamL","ShamR","Hour4L","Hour4R",
                "Hour12L","Hour12R","Day2L","Day2R",
                "Day14L","Day14R","Week6L","Week6R"),
  condition = c("Sham","Sham","Hour4","Hour4",
                "Hour12","Hour12","Day2","Day2",
                "Day14","Day14","Week6","Week6"),
  side      = c("L","R","L","R","L","R","L","R","L","R","L","R"),
  chip      = c("0005295","0005295","0005161","0005161",
                "0005295","0005295","0005161","0005161",
                "0005161","0005161","0005295","0005295"),
  inner_dir = c(
    "ShamL-output-XETG00063__0005295__Region_5__20230717__191520",
    "ShamR-output-XETG00063__0005295__Region_6__20230717__191520",
    "Hour4L-output-XETG00063__0005161__Region_5__20230717__191520",
    "Hour4R-output-XETG00063__0005161__Region_6__20230717__191520",
    "Hour12L-output-XETG00063__0005295__Region_3__20230717__191520",
    "Hour12R-output-XETG00063__0005295__Region_4__20230717__191520",
    "Day2L-output-XETG00063__0005161__Region_3__20230717__191519",
    "Day2R-output-XETG00063__0005161__Region_4__20230717__191519",
    "Day14L-output-XETG00063__0005161__Region_1__20230717__191519",
    "Day14R-output-XETG00063__0005161__Region_2__20230717__191519",
    "Week6L-output-XETG00063__0005295__Region_1__20230717__191520",
    "Week6R-output-XETG00063__0005295__Region_2__20230717__191520"
  ),
  stringsAsFactors = FALSE
)

# ── Step 1: Load each Xenium sample ───────────────────────────────────────────
message("==> Loading Xenium samples ...")

load_xenium_sample <- function(row) {
  data_dir <- file.path(XEN_BASE, row$inner_dir)
  message(sprintf("  Loading %s ...", row$sample_id))

  obj <- LoadXenium(data_dir, fov = "fov", segmentations = "cell")

  # QC filter matching published IRI pipeline
  obj <- subset(obj, subset = nCount_Xenium > 10 & nFeature_Xenium > 3)

  obj$sample_id <- row$sample_id
  obj$condition <- row$condition
  obj$side      <- row$side
  obj$chip      <- row$chip
  obj$gsm       <- row$gsm

  log_mem(paste("after loading", row$sample_id))
  obj
}

sample_list <- vector("list", nrow(sample_meta))
for (i in seq_len(nrow(sample_meta))) {
  sample_list[[i]] <- load_xenium_sample(sample_meta[i, ])
}
names(sample_list) <- sample_meta$sample_id

cell_counts <- sapply(sample_list, ncol)
message("  Cell counts per sample:")
print(cell_counts)
message(sprintf("  Total cells: %d", sum(cell_counts)))

# ── Step 2: Merge all 12 samples ──────────────────────────────────────────────
message("==> Merging all 12 samples ...")
object <- merge(
  x            = sample_list[[1]],
  y            = sample_list[-1],
  add.cell.ids = names(sample_list),
  project      = "Xenium_UUO"
)
rm(sample_list); gc()
log_mem("after merge")

DefaultAssay(object) <- "Xenium"

# ── Step 3: Normalize ─────────────────────────────────────────────────────────
# Log-normalization with scale.factor=100, matching published IRI pipeline
message("==> Normalizing ...")
object <- NormalizeData(object,
  normalization.method = "LogNormalize",
  scale.factor         = 100
)
object <- JoinLayers(object)
log_mem("after normalization")

# ── Step 4: PCA on all cells ──────────────────────────────────────────────────
# Use all Xenium genes — panel is pre-curated, no HVG selection needed
# No sketching — ~300 genes makes PCA fast even on millions of cells
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

# ── Step 5: Harmony batch correction ──────────────────────────────────────────
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

# ── Step 6: Cluster and UMAP on all cells ─────────────────────────────────────
message("==> Clustering and UMAP ...")
object <- FindNeighbors(object, reduction = "harmony", dims = 1:30)
object <- FindClusters(object, resolution = 0.3, cluster.name = "seurat_clusters")
object <- RunUMAP(object,
  reduction      = "harmony",
  dims           = 1:30,
  reduction.name = "umap"
)
log_mem("after clustering")
# save the object after clustering and UMAP, before label transfer
saveRDS(object, file.path(OUT_DIR, "Xenium_harmony_integrated_prelabeltransfer.rds"))
message("  Saved: ", file.path(OUT_DIR, "Xenium_harmony_integrated_prelabeltransfer.rds"))


# helper: prepare a reference for label transfer
# - subsets to shared genes with Xenium panel
# - scales and runs PCA on the data layer
prepare_reference <- function(ref_path, label) {
  message(sprintf("==> Loading %s reference ...", label))
  ref <- readRDS(ref_path)
  DefaultAssay(ref) <- "originalexp"

  # Reference rownames are Ensembl IDs (e.g. "ENSMUSG00000026315"); the
  # actual gene symbol is embedded in meta.features$feature_name, as either
  # "Symbol" or "Symbol_ENSMUSG..." (suffixed only for duplicated symbols).
  # Swap rownames to symbols so they can be matched against the Xenium panel.
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

# ── Step 7: Label transfer from MKA reference ─────────────────────────────────
mka      <- prepare_reference(REF_MKA, "MKA")
ref_mka  <- mka$ref
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

# ── Step 8: Label transfer from Lake snRNA-seq reference ──────────────────────
lake       <- prepare_reference(REF_LAKE, "Lake")
ref_lake   <- lake$ref
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

# ── Step 9: Visualize and save ────────────────────────────────────────────────
message("==> Saving plots ...")
pdf(OUT_PDF, width = 16, height = 7)

p_cluster <- DimPlot(object,
  reduction = "umap", group.by = "seurat_clusters",
  label = TRUE, repel = TRUE, pt.size = 0.1
) + ggtitle("Harmony clusters")

p_mka <- DimPlot(object,
  reduction = "umap", group.by = "predicted.MKA_celltype",
  label = TRUE, repel = TRUE, pt.size = 0.1
) + ggtitle("MKA: author_cell_type")

p_lake1 <- DimPlot(object,
  reduction = "umap", group.by = "predicted.Lake_SubclassLevel1",
  label = TRUE, repel = TRUE, pt.size = 0.1
) + ggtitle("Lake: SubclassLevel1 (broad)")

p_lake3 <- DimPlot(object,
  reduction = "umap", group.by = "predicted.Lake_SubclassLevel3",
  label = TRUE, repel = TRUE, pt.size = 0.1
) + ggtitle("Lake: SubclassLevel3 (detailed)")

print(p_cluster | p_mka)
print(p_lake1 | p_lake3)

p_sample <- DimPlot(object,
  reduction = "umap", group.by = "sample_id",
  label = FALSE, pt.size = 0.1
) + ggtitle("Sample identity")

p_cond <- DimPlot(object,
  reduction = "umap", group.by = "condition",
  label = FALSE, pt.size = 0.1
) + ggtitle("Condition (UUO time-point)")

p_side <- DimPlot(object,
  reduction = "umap", group.by = "side",
  label = FALSE, pt.size = 0.1
) + ggtitle("Side (L=obstructed, R=contralateral)")

print(p_sample | p_cond | p_side)

# Spatial plots per sample
for (i in seq_len(nrow(sample_meta))) {
  samp    <- sample_meta$sample_id[i]
  cells_s <- WhichCells(object, expression = sample_id == samp)
  sub_obj <- subset(object, cells = cells_s)
  tryCatch({
    print(
      ImageDimPlot(sub_obj, fov = "fov", group.by = "predicted.MKA_celltype",
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

message("==> Done.")
log_mem("final")
