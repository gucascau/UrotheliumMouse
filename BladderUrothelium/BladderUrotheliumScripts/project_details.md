# Bladder Urothelium scRNA-seq Integration Project

**Species:** Mouse (Mus musculus)  
**Last updated:** 2026-05-05  
**Working directory:** `/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/BladderUrothelium/`

---

## Datasets

| Sample ID | GSM Accession | Accession Series | Technology | Condition | Reference |
|-----------|---------------|------------------|-----------|-----------|-----------|
| BladderNormal1 | GSM3723360 | GSE129845 | 10X v2 | Normal bladder | Yu et al. 2019 (PMID 31462402) |
| BladderNormal2 | GSM3723361 | GSE129845 | 10X v2 | Normal bladder | Yu et al. 2019 (PMID 31462402) |
| BladderUro1 | GSM4970435 | GSE163029 | 10X | Healthy urothelium | Liu et al. 2021 (PMID 33538002) |
| BladderWT1 | GSM4201633 | GSE141348 | 10X | WT mouse bladder | — |
| BladderB8W1 | GSM5014059 | GSE164557 | 10X | Healthy 8 weeks | — |
| BladderB8W2 | GSM5014060 | GSE164557 | 10X | Healthy 8 weeks | — |
| BladderH48_1 | GSM5014061 | GSE164557 | 10X | Acute CPP injury (48h) | — |
| BladderH48_2 | GSM5014062 | GSE164557 | 10X | Acute CPP injury (48h) | — |
| BladderD11_1 | GSM5014063 | GSE164557 | 10X | Chronic CPP injury (11 days) | — |
| BladderD11_2 | GSM5014064 | GSE164557 | 10X | Chronic CPP injury (11 days) | — |

**Raw data formats:**
- GSE129845 — standard 10X MTX (gzipped barcodes/genes/matrix)
- GSE163029 — standard 10X MTX (gzipped barcodes/genes/matrix)
- GSE141348 (GSM4201633) — 10X v3 barnyard features (mm10 + GRCh38); mouse genes extracted by filtering `mm10___ENSMUSG` features
- GSE164557 — dense tab-delimited XLS (Gene_ID | Symbol | cells), converted to sparse matrix

---

## Cell Counts

### After QC + Doublet Removal (per-sample h5ad)

| Sample ID | Condition | Cells (post-QC) | Genes |
|-----------|-----------|-----------------|-------|
| BladderUro1 | Healthy urothelium | 17,273 | 17,119 |
| BladderH48_2 | Acute CPP 48h | 13,453 | 18,056 |
| BladderD11_1 | Chronic CPP 11d | 10,363 | 16,316 |
| BladderB8W2 | Healthy 8wk | 9,063 | 16,499 |
| BladderB8W1 | Healthy 8wk | 7,872 | 15,213 |
| BladderWT1 | WT bladder | 6,922 | 18,597 |
| BladderNormal2 | Normal bladder | 6,799 | 17,003 |
| BladderNormal1 | Normal bladder | 6,753 | 16,744 |
| BladderH48_1 | Acute CPP 48h | 6,570 | 15,065 |
| BladderD11_2 | Chronic CPP 11d | 5,162 | 15,889 |
| **TOTAL** | | **90,230** | 20,707 (shared) |

### After scVI Integration (all cells)

| Condition | Cells |
|-----------|-------|
| Acute CPP 48h | 20,023 |
| Healthy urothelium | 17,273 |
| Healthy 8wk | 16,935 |
| Chronic CPP 11d | 15,525 |
| Normal bladder | 13,552 |
| WT bladder | 6,922 |
| **TOTAL** | **90,230** |

### After Urothelial Cell Extraction

| Sample ID | Condition | Urothelial Cells |
|-----------|-----------|-----------------|
| BladderUro1 | Healthy urothelium | 17,273 |
| BladderH48_2 | Acute CPP 48h | 12,848 |
| BladderD11_1 | Chronic CPP 11d | 10,358 |
| BladderB8W2 | Healthy 8wk | 8,972 |
| BladderB8W1 | Healthy 8wk | 7,331 |
| BladderH48_1 | Acute CPP 48h | 6,423 |
| BladderNormal2 | Normal bladder | 5,549 |
| BladderNormal1 | Normal bladder | 5,506 |
| BladderD11_2 | Chronic CPP 11d | 5,153 |
| BladderWT1 | WT bladder | 4,302 |
| **TOTAL** | | **83,715** |

| Condition | Urothelial Cells |
|-----------|-----------------|
| Acute CPP 48h | 19,271 |
| Healthy urothelium | 17,273 |
| Healthy 8wk | 16,303 |
| Chronic CPP 11d | 15,511 |
| Normal bladder | 11,055 |
| WT bladder | 4,302 |

---

## Processing Pipeline

| Step | Script | Description | Output |
|------|--------|-------------|--------|
| 01 | `01_load_bladder.R` | Load raw data → individual Seurat objects | `seurat_objects/*_seurat.rds` |
| 02 | `02_qc_bladder.R` | QC filtering (nFeature 200–8000, nCount ≥300, pct_mt ≤40%) + DoubletFinder | `seurat_objects/*_qc.rds` |
| 03 | `03_export_h5ad_bladder.R` | Export QC'd Seurat → h5ad | `qc_h5ad/*.h5ad` |
| 04 | `04_harmony_bladder.R` | Seurat+Harmony integration (all cells; batch = sample_id + technology) | `integration_output/bladder_harmony_integrated.rds` |
| 05b | `05b_scvi_fullgene.py` | Full-gene scVI integration (3,000 HVGs for training; full ~20k gene output); n_latent=20, n_layers=2, n_hidden=256, 400 epochs | `output/BladderUrothelium_allcells_scvi.h5ad` |
| 06 | `06_extract_uro.py` | Urothelial cell extraction (three-arm OR gate: umbrella UPK, basal ≥2/3 Krt5/Krt14/Trp63, intermediate Epcam+Krt7+Krt8+Foxa1/Gata3) | `output/BladderUrothelium_uro_cells_scvi.h5ad` |
| 07 | `07_convert_to_rds.R` | Convert urothelial h5ad → Seurat RDS | `output/BladderUrothelium_uro_cells_scvi.rds` |
| 08 | `08_harmony_recluster.R` | Harmony reclustering of urothelial cells only (20 PCs, res=0.5) | `output/BladderUrothelium_uro_cells_harmony_integrated.rds` |

### QC Thresholds (Script 02)
- `min_features` = 200, `max_features` = 8,000
- `min_counts` = 300
- `max_pct_mt` = 40%
- DoubletFinder rate: `min(0.008 × n_cells/1000, 0.25)`

### scVI Parameters (Script 05b)
- n_HVG = 3,000 (batch-aware selection)
- n_latent = 20, n_layers = 2, n_hidden = 256
- dropout = 0.1, dispersion = gene-batch, likelihood = NB
- 400 epochs, batch_size = 256, lr = 1e-3

---

## Output Files

| File | Location | Description |
|------|----------|-------------|
| `*_seurat.rds` | `seurat_objects/` | Raw per-sample Seurat objects (10 files) |
| `*_qc.rds` | `seurat_objects/` | QC + doublet-filtered Seurat objects (10 files) |
| `*.h5ad` | `qc_h5ad/` | QC'd per-sample AnnData objects (10 files) |
| `bladder_harmony_integrated.rds` | `integration_output/` | Seurat Harmony integration (all cells) |
| `bladder_scvi_integrated.h5ad` | `integration_output/` | scVI integration (HVG subset, all cells) |
| `BladderUrothelium_allcells_scvi.h5ad` | `output/` | Full-gene scVI integration (90,230 cells × 20,707 genes) |
| `BladderUrothelium_uro_cells_scvi.h5ad` | `output/` | Urothelial cells only (83,715 cells × 20,707 genes) |
| `BladderUrothelium_uro_cells_scvi.rds` | `output/` | Urothelial cells as Seurat RDS |
| `bladder_fullgene_scvi_model/` | `output/` | Saved scVI model directory |

---

## Key Notes

- **GSM4201633 (BladderWT1)** was aligned to a combined mm10+GRCh38 barnyard reference; only `mm10___ENSMUSG` features were retained, and the `mm10___` prefix was stripped to match gene symbol conventions across samples.
- **GSE164557 samples** (B8W, H48, D11) use dense tab-delimited XLS format; cells were converted to sparse matrix on load.
- **Urothelial extraction** used a three-arm OR gate (umbrella UPK, basal ≥2/3 keratin, intermediate Epcam+TF). This recovered 83,715 / 90,230 cells (92.8%), indicating the dataset is predominantly urothelial-enriched.
- **scVI trained on HVG subset** (3,000 genes) but latent embeddings were applied to the full-gene AnnData for downstream analysis.
- SLURM job `7931283` ran the full pipeline (steps 05b → 08) on 2026-04-27/28.
