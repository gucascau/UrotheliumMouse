---
name: UsedSingleCells dataset inventory
description: All scRNA-seq datasets in the UsedSingleCells working directory — accessions, formats, cell counts, conditions
type: project
originSessionId: 9c109388-02fc-4a35-a3e8-50c8d17651ce
---
Working directory: /vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells

## Dataset Summary

| Sample(s) | Accession | Technology | Format | Species | Condition | Has annotation | ~Cells |
|-----------|-----------|-----------|--------|---------|-----------|---------------|--------|
| Healthy, UUO | GSE119531 | Drop-seq | Dense DGE txt + cell annotation TSV | mouse | Healthy / UUO | YES (CellType col) | 11,395 / 6,147 |
| kid | GSE108291 | 10X v2 | MTX (unfiltered, 2.2M barcodes) | mouse | Normal kidney | No | ~few k after filter |
| org, org4 | GSE108291 | 10X v2 | MTX | **human** | Kidney organoid | No | ~1,421 (org4) |
| sham, uuo2, uuo7, ruuo | GSM4151577–4580 | 10X v2 | MTX (filtered) | mouse | Sham/UUO-day2/day7/rUUO (Kirita 2020) | No | 5938/2739/4504/3695 |
| CS1KKY, CU1KKY, CU2KKY | GSM5333084–86 | 10X v3 (CR4) | MTX (unfiltered, 6.8M barcodes) | mouse | CtrlSham / ChronicUUO | No | ~few k after filter |
| Sample4, Sample5 | GSM3723360–61 | 10X v2 | MTX (filtered) | mouse | Unknown | No | 7820 / 7975 |
| sample3–6 | GSM8417119–22 | 10X Multiome | MTX with RNA+ATAC (181k features: 32,285 RNA + 148,972 ATAC) | mouse | Unknown | No | 8125–20000 |
| sci_kidney | GSE119945 | sci-RNA-seq | MatrixMarket + cell_annotate.csv + gene_annotate.csv | mouse | Kidney development | YES (Main_Cluster, Sub_Cluster) | ~2M raw, filter to non-doublets |
| DU1, DU2 | GSE264184 | Unknown | Dense gene×cell TSV (pre-normalized float values) | mouse | Disease (DU) | No | ~14,559 total |
| KudoUUO | RDS | 10X | Seurat RDS (already Harmony-integrated) | mouse | UUO comparison | YES | unknown |
| MKA | RDS | mixed | SCE (zellkonverter) – Mouse Kidney Atlas | mouse | Atlas | YES | unknown |
| ChenSpatial | RDS | spatial | SCE (zellkonverter) – Chen 2025 NatGenet | mouse | Sex-specific lifespan | YES | unknown |
| LakesnRNA | RDS | snRNA-seq | SCE (zellkonverter) – Lake 2025 bioRxiv | mouse | Reference | YES | unknown |

## BladderUrothelium Datasets (BladderUrothelium/Raw/)

| Sample(s) | Accession | Technology | Format | Species | Condition | Has annotation | ~Cells |
|-----------|-----------|-----------|--------|---------|-----------|---------------|--------|
| GSM4970435 | GSE163029 | 10X Genomics | MTX (barcodes/genes/matrix gz) | mouse | Healthy urothelium | No | 18,917 |
| GSM4201633 | GSE141348 | 10X Genomics | MTX (barcodes/features/matrix gz) | mouse | WT mouse bladder | No | ~? |
| B8W-1, B8W-2 | GSE164557 | 10X Genomics | Dense gene×cell XLS gz | mouse | Healthy bladder (8wk) | No | ~? each |
| H48-1, H48-2 | GSE164557 | 10X Genomics | Dense gene×cell XLS gz | mouse | Acute CPP injury (48h) | No | ~? each |
| D11-1, D11-2 | GSE164557 | 10X Genomics | Dense gene×cell XLS gz | mouse | Chronic CPP injury (11 days) | No | ~? each |
| Sample4, Sample5 | GSE129845 / GSM3723360-61 | 10X v2 | MTX (symlinks to parent dir) | mouse | Normal bladder homogenate | No | 7820 / 7975 |

Paper refs:
- GSE163029 (PMID 33538002): Liu et al. 2021 — mouse bladder urothelium 8 subpopulations
- GSE164557: Bladder urothelium & stroma heterogeneity in acute and chronic cyclophosphamide injury (10X, NovaSeq)
- GSE129845 (PMID 31462402): Yu et al. 2019 — human & mouse bladder map (mouse samples = Sample4/5)

## Kidney Datasets (UsedSingleCells/)

| Sample(s) | Accession | Technology | Format | Species | Condition | Has annotation | ~Cells |
|-----------|-----------|-----------|--------|---------|-----------|---------------|--------|
| Healthy, UUO | GSE119531 | Drop-seq | Dense DGE txt + cell annotation TSV | mouse | Healthy / UUO | YES (CellType col) | 11,395 / 6,147 |
| kid | GSE108291 | 10X v2 | MTX (unfiltered, 2.2M barcodes) | mouse | Normal kidney | No | ~few k after filter |
| org, org4 | GSE108291 | 10X v2 | MTX | **human** | Kidney organoid | No | ~1,421 (org4) |
| sham, uuo2, uuo7, ruuo | GSM4151577–4580 | 10X v2 | MTX (filtered) | mouse | Sham/UUO-day2/day7/rUUO (Kirita 2020) | No | 5938/2739/4504/3695 |
| CS1KKY, CU1KKY, CU2KKY | GSM5333084–86 | 10X v3 (CR4) | MTX (unfiltered, 6.8M barcodes) | mouse | CtrlSham / ChronicUUO | No | ~few k after filter |
| sample3–6 | GSM8417119–22 | 10X Multiome | MTX with RNA+ATAC (181k features: 32,285 RNA + 148,972 ATAC) | mouse | Unknown | No | 8125–20000 |
| sci_kidney | GSE119945 | sci-RNA-seq | MatrixMarket + cell_annotate.csv + gene_annotate.csv | mouse | Kidney development | YES (Main_Cluster, Sub_Cluster) | ~2M raw, filter to non-doublets |
| DU1, DU2 | GSE264184 | Unknown | Dense gene×cell TSV (pre-normalized float values) | mouse | Disease (DU) | No | ~14,559 total |
| KudoUUO | RDS | 10X | Seurat RDS (already Harmony-integrated) | mouse | UUO comparison | YES | unknown |
| MKA | RDS | mixed | SCE (zellkonverter) – Mouse Kidney Atlas | mouse | Atlas | YES | unknown |
| ChenSpatial | RDS | spatial | SCE (zellkonverter) – Chen 2025 NatGenet | mouse | Sex-specific lifespan | YES | unknown |
| LakesnRNA | RDS | snRNA-seq | SCE (zellkonverter) – Lake 2025 bioRxiv | mouse | Reference | YES | unknown |

## Key Notes
- GSE108291 org/org4 use HUMAN Ensembl IDs (ENSG) — need ortholog mapping to integrate with mouse data
- GSM5333084–86 and GSE108291 kid are UNFILTERED — must run emptyDrops before creating Seurat objects
- GSM8417119–22 are multiome — filter features to "Gene Expression" only for RNA integration
- GSE264184 contains non-integer (normalized) values — skip NormalizeData, load directly as normalized data
- GSE119945 sci-RNA-seq has 2M barcodes — filter to annotated non-doublet cells before loading
