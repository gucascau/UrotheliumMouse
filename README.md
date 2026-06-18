# Mouse Urothelium Single-Cell RNA-seq Integration Pipeline

Single-cell and single-nucleus RNA-seq integration pipeline for mouse urothelium datasets across kidney, bladder, and ureter tissues. Integrates 26 datasets from multiple platforms and tissue types using Harmony and scVI/scANVI.

## Datasets

### Kidney

| Sample ID | Accession | Technology | Condition |
|---|---|---|---|
| KidneyHealthy1 | GSE119531 | Drop-seq | Healthy |
| KidneyUUO1 | GSE119531 | snRNA-seq | UUO |
| KidneyE18_5_1 | GSE108291 | 10X v2 | Normal (E18.5) |
| KidneyHealthy2 | GSM4151577 | 10X v2 | Sham |
| KidneyUUO2 | GSM4151578 | 10X v2 | UUO day 2 |
| KidneyUUO3 | GSM4151579 | 10X v2 | UUO day 7 |
| KidneyrUUO1 | GSM4151580 | 10X v2 | rUUO |
| KidneyHealthy3 | GSM5333084 | 10X v3.1 | Sham |
| KidneyUUO4 | GSM5333085 | 10X v3.1 | UUO 10 days |
| KidneyUUO5 | GSM5333086 | 10X v3.1 | UUO 10 days |
| KidneyHealthy4 | GSM8417119 | snRNA-seq | Sham |
| KidneyHealthy5 | GSM8417120 | snRNA-seq | Sham |
| KidneyUUO6 | GSM8417121 | snRNA-seq | UUO 4 days |
| KidneyTET2UUO | GSM8417122 | snRNA-seq | TET2KO UUO 4 days |
| EmbryosE9_5ToE13_5 | GSE119945 | sci-RNA-seq3 | Development |
| KidneyUUO7 | GSE264184 | Dense matrix | UUO 10 days |
| KidneyUUO8 | GSE264184 | Dense matrix | UUO 10 days |
| MKA | MKARDS | Mixed | Atlas |
| ChenSpatial | ChenRDS | Mixed | Sex-specific |
| LakesnRNA | LakesnRDS | Mixed | Reference |

### Bladder

| Sample ID | Accession | Technology | Condition | Processed File |
|---|---|---|---|---|
| BladderUrothelium | GSE129845 / GSE163029 / GSE164557 / GSM4201633 | 10X | Healthy bladder | `BladderUrothelium_uro_cells_scvi.rds` |
| BladderHomogenate1 | In-house | 10X | Bladder organoid | `BladderHomogenate1_qc.rds` |
| BladderHomogenate2 | In-house | 10X | Bladder organoid | `BladderHomogenate2_qc.rds` |

### Ureter

| Sample ID | Accession | Technology | Condition | Processed File |
|---|---|---|---|---|
| UreterOrganoid | In-house | 10X | Ureter organoid | `MouseUreterRecon1_qc.rds` |

### Organoids

| Sample ID | Accession | Technology | Condition | Processed File |
|---|---|---|---|---|
| KudoUUOUrothelium | In-house | 10X | UUO urothelium organoid | `KudoUUOUrothelium_qc.rds` |

## Pipeline Overview

```
01_load_datasets.R          # Load raw data → Seurat objects
02_qc_and_filter.R          # QC, filtering, DoubletFinder
03_integrate_harmony.R      # Merge, normalize, HVG selection → saves checkpoint
04_integrate_harmony.R      # ScaleData, PCA, Harmony, UMAP (runs from checkpoint)
03b_export_for_scvi.R       # Export raw counts → scvi_input.h5ad
04_scvi_integration.py      # scVI / scANVI integration (Python)
```

## Requirements

### R packages
- Seurat v5, harmony, DoubletFinder
- SingleCellExperiment, DropletUtils
- org.Mm.eg.db (Ensembl → gene symbol mapping)
- ggplot2, dplyr, patchwork, Matrix

### Python packages (conda env: `cell2loc_env`)
- scvi-tools >= 1.2.1
- scanpy >= 1.9
- torch 2.5.1+cu124
- leidenalg, igraph, anndata

## Usage (HPC / SLURM)

**Step 1 — Load datasets**
```bash
sbatch submit_01_load.sh
```

**Step 2 — QC and doublet removal** (array job, one task per sample)
```bash
sbatch submit_02_qc_array.sh
```

**Step 3 — Merge, normalize, HVG selection**
```bash
sbatch submit_03_integrate.sh
# Outputs: integration_output/merged_normalized.rds
#          integration_output/hvg_list.rds
```

**Step 4a — Harmony integration** (run after Step 3)
```bash
sbatch submit_04_harmony.sh
# Output: integration_output/merged_harmony_integrated.rds
```

**Step 4b — Export for scVI** (run after Step 3, can run in parallel with 4a)
```bash
sbatch submit_03b_export.sh
# Output: integration_output/scvi_input.h5ad
```

**Step 4c — scVI / scANVI integration** (run after Step 4b)
```bash
sbatch submit_04_scvi.sh
# Outputs: integration_output/scvi_integrated.h5ad
#          integration_output/scvi_model/
#          integration_output/scanvi_model/
```

## Key Design Decisions

- **Ensembl ID conversion**: ChenSpatial, LakesnRNA, MKA, and KudoUUOUrothelium datasets contain Ensembl IDs; these are converted to gene symbols using `org.Mm.eg.db` before merging. Duplicates are resolved by keeping the highest-expressing gene.
- **Batch correction**: Harmony corrects for both `sample_id` and `technology` (Drop-seq, 10X v2/v3, snRNA-seq, sci-RNA-seq3).
- **scVI exclusions**: ChenSpatial, LakesnRNA, MKA, KudoUUOUrothelium, KidneyUUO7/8 are excluded from scVI because they lack raw integer counts (pre-normalized matrices).
- **Memory**: Full dataset (~2M+ cells) requires 512 GB RAM (himem partition). scVI training uses GPU.
