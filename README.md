# Mouse Urothelium Single-Cell RNA-seq Integration Pipeline

Single-cell and single-nucleus RNA-seq integration pipeline for mouse urothelium datasets across kidney, bladder, and ureter tissues. Integrates 34 datasets from multiple platforms and tissue types using Harmony and scVI/scANVI.

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

10 scRNA-seq samples integrated into `BladderUrothelium_uro_cells_scvi.rds` (Yu et al. 2019 PMID:31462402; Liu et al. 2021 PMID:33538002; GSE164557; GSM4201633).

| Sample ID | Accession | Technology | Condition |
|---|---|---|---|
| BladderNormal1 | GSM3723360 | 10X v2 | Normal bladder |
| BladderNormal2 | GSM3723361 | 10X v2 | Normal bladder |
| BladderUro1 | GSM4970435 | 10X | Healthy urothelium |
| BladderWT1 | GSM4201633 | 10X | WT mouse bladder |
| BladderB8W1 | GSM5014059 | 10X | Healthy 8 wks |
| BladderB8W2 | GSM5014060 | 10X | Healthy 8 wks |
| BladderH48_1 | GSM5014061 | 10X | Acute CPP 48h |
| BladderH48_2 | GSM5014062 | 10X | Acute CPP 48h |
| BladderD11_1 | GSM5014063 | 10X | Chronic CPP 11 days |
| BladderD11_2 | GSM5014064 | 10X | Chronic CPP 11 days |
| BladderHomogenate1 | GSM3827175 | 10X | Bladder organoid |
| BladderHomogenate2 | GSM3827176 | 10X | Bladder organoid |

### Ureter

| Sample ID | Accession | Technology | Condition | Processed File |
|---|---|---|---|---|
| UreterOrganoid | In-house | 10X | Ureter organoid | `MouseUreterRecon1_qc.rds` |

### Organoids

| Sample ID | Accession | Technology | Condition | Processed File |
|---|---|---|---|---|
| KudoUUOUrothelium | GSE190887 / In-house | PIPseq / sci-RNA-seq3 | UUO urothelium organoid + in vivo | `KudoUUOUrothelium_qc.rds` |

## Pipeline Overview

```text
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

### Step 1 — Load datasets

```bash
sbatch submit_01_load.sh
```

### Step 2 — QC and doublet removal (array job, one task per sample)

```bash
sbatch submit_02_qc_array.sh
```

### Step 3 — Merge, normalize, HVG selection

```bash
sbatch submit_03_integrate.sh
# Outputs: integration_output/merged_normalized.rds
#          integration_output/hvg_list.rds
```

### Step 4a — Harmony integration (run after Step 3)

```bash
sbatch submit_04_harmony.sh
# Output: integration_output/merged_harmony_integrated.rds
```

### Step 4b — Export for scVI (run after Step 3, can run in parallel with 4a)

```bash
sbatch submit_03b_export.sh
# Output: integration_output/scvi_input.h5ad
```

### Step 4c — scVI / scANVI integration (run after Step 4b)

```bash
sbatch submit_04_scvi.sh
# Outputs: integration_output/scvi_integrated.h5ad
#          integration_output/scvi_model/
#          integration_output/scanvi_model/
```

## Urothelium Extraction Criteria (4-arm gate)

Applied in `UrotheliumIntegrationScripts/02b_plot_AllUrothelium_markers_AllscRNA.R` after initial integration.
A cell is **kept** if it passes **any arm** AND the **kidney score filter**.

### Arms

| Arm | Markers | Logic | Tissues |
|-----|---------|-------|---------|
| **Umbrella** | Upk1a, Upk1b, Upk2, Upk3a, Upk3b | any > 0 | All |
| **Basal** | Krt5, Krt14, Trp63 | ≥ 2 of 3 > 0 | All |
| **Intermediate** | Epcam, Krt7, Krt8, Foxa1/Gata3 | all four expressed simultaneously | All |
| **Pan-keratin** | Krt8, Krt18, Krt19 | any > 0 | **Non-kidney only** |

> **Why pan-keratin is excluded from kidney:**
> Krt18 is expressed in kidney proximal tubular cells, so it cannot safely gate
> kidney urothelium on its own. Bladder and ureter have no tubular contamination,
> so Krt8/Krt18/Krt19 are safe universal markers there. Ureter urothelium also
> expresses lower uroplakins than bladder, and the strict intermediate arm
> (requiring all 4 markers simultaneously) loses genuine cells to scRNA-seq
> dropout — the pan-keratin arm recovers them.

### Kidney-specific filter (KUDO qc cells only)

```text
KidneyEpiScore = AddModuleScore(tubular markers)
Keep if KidneyEpiScore1 < 0.2
```

Tubular markers scored:

| Cell type | Markers |
|-----------|---------|
| Proximal tubule | Slc34a1, Lrp2, Cubn |
| TAL | Umod, Slc12a1 |
| DCT | Slc12a3 |
| Collecting duct principal | Aqp2, Aqp3, Scnn1g |
| Intercalated cells | Atp6v1b1, Slc4a1, Foxi1 |

> **Why only KUDO qc cells?**
> scVI-processed kidney in vivo cells were already filtered by the upstream
> `03_extract_uro.py` pipeline. Applying AddModuleScore to a mixed-tissue
> dataset shifts the background and inflates scores for genuine urothelium.

### Final gate logic

```text
keep = (umbrella | basal | intermediate | pan_krt) & kidney_score_ok
```

## Key Design Decisions

- **Ensembl ID conversion**: ChenSpatial, LakesnRNA, MKA, and KudoUUOUrothelium datasets contain Ensembl IDs; these are converted to gene symbols using `org.Mm.eg.db` before merging. Duplicates are resolved by keeping the highest-expressing gene.
- **Batch correction**: Harmony corrects for both `sample_id` and `technology` (Drop-seq, 10X v2/v3, snRNA-seq, sci-RNA-seq3).
- **scVI exclusions**: ChenSpatial, LakesnRNA, MKA, KudoUUOUrothelium, KidneyUUO7/8 are excluded from scVI because they lack raw integer counts (pre-normalized matrices).
- **Memory**: Full dataset (~2M+ cells) requires 512 GB RAM (himem partition). scVI training uses GPU.
