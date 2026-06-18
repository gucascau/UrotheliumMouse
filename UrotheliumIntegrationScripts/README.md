# Mouse Urothelium Atlas — Analysis Pipeline

**Project:** Multi-tissue, multi-chemistry integration of mouse urothelium cells  
**Tissues:** Kidney · Bladder · Ureter (in vivo + organoid)  
**Contact:** xin.wang@nationwidechildrens.org

---

## Data Sources

Six input RDS files are merged in `01_integrate`:

| File | Tissue | Context | Chemistry | Normalisation |
|------|--------|---------|-----------|---------------|
| `BladderUrothelium_uro_cells_scvi.rds` | Bladder | In vivo | 10X | scVI log-norm (in counts layer) |
| `RenalUrothelium_uro_cells_fullgene_scvi.rds` | Kidney | In vivo | 10X | scVI log-norm (in counts layer) |
| `BladderHomogenate1_qc.rds` (GSM3827175) | Bladder | Organoid | 10X | Raw counts → NormalizeData |
| `BladderHomogenate2_qc.rds` (GSM3827176) | Bladder | Organoid | 10X | Raw counts → NormalizeData |
| `KudoUUOUrothelium_qc.rds` | Kidney | Organoid + sci-RNA-seq3 | PIPseq / sci-RNA-seq3 | Raw counts → NormalizeData |
| `MouseUreterRecon1_qc.rds` | Ureter | Organoid | 10X | Raw counts → NormalizeData |

The KUDO object contains two sub-populations:
- **PIPseq organoid** (GOF / LOF / Vehicle / Yoda) — kidney urothelium organoid
- **sci-RNA-seq3 in vivo** (Health / UUO_Day2 through Day14) — kidney in vivo (GSE190887)

---

## Pipeline Overview

```
01_integrate_AllUrothelium_harmony.R
│   Merge 6 sources → 3 000 HVGs → ScaleData(pct_mt) → PCA(20) →
│   Harmony(sample_id + technology + source) → Clusters(r=0.5) → UMAP
│   OUTPUT: AllUrothelium_harmony_integrated.rds
│
└─► 02b_plot_AllUrothelium_markers_AllscRNA.R
    │   Metadata fixes → 4-arm urothelium gate → re-integrate gated cells:
    │   3 000 HVGs → ScaleData → PCA(50) → Harmony(sample_id + technology) →
    │   Clusters(r=0.5) → UMAP → marker plots
    │   OUTPUT: AllUrothelium_markers_gated.rds
    │           AllUrothelium_gated_UMAP.rds          ← main downstream input
    │
    ├─► 03b_metadata_correlations_AllscRNAUrothelium.R
    │       Drop snRNA-seq → avg-expression correlations (3 gene modes) →
    │       Cramer's V metadata associations → cluster-composition correlations
    │       OUTPUT: AllUrothelium_gated_UMAP_nosnRNA.rds
    │                   │
    │                   └─► 02c_UMAP_plots_AllscRNA.R
    │                           Re-cluster → UMAP → FeaturePlot / DotPlot / DimPlot
    │                           OUTPUT: AllUrothelium_nosnRNA_UMAP.rds
    │
    ├─► 04a_export_for_scvi.R
    │       Batch-aware 4 000 HVGs → export MTX + metadata
    │       OUTPUT: AllUrothelium_scvi_input/
    │                   │
    │                   └─► 04b_scvi_AllUrothelium.py
    │                           scVI(n_latent=30, n_layers=2, normal likelihood,
    │                           max_epochs=400) → kNN → UMAP → Leiden clusters
    │                           OUTPUT: AllUrothelium_scvi.h5ad
    │
    ├─► 06_compare_urothelium_tissues.R
    │       DE (FindMarkers/Wilcoxon) + volcano + GO/KEGG + PROGENy
    │       5 comparisons: kidney vs bladder (all), organoid pairwise,
    │       kidney invivo vs organoid, healthy-only variants
    │       OUTPUT: AllUrothelium_tissue_comparison/
    │
    └─► 07_shared_function_urothelium.R
            One-vs-rest FindMarkers per Category → UpSet plot →
            core genes (upregulated in ≥ 3 groups) → GO/KEGG →
            pan-urothelial module score on UMAP
            OUTPUT: AllUrothelium_shared_function/
```

---

## Urothelium Extraction Criteria (4-arm gate)

Applied in `02b_plot_AllUrothelium_markers_AllscRNA.R` after initial integration.  
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

```
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

```
keep = (umbrella | basal | intermediate | pan_krt) & kidney_score_ok
```

---

## Post-gate Metadata Fixes

Applied in `02b` before and after gating:

| Issue | Fix |
|-------|-----|
| `tissue = NA` for most cells | Fill with "kidney"; GSM3827175/76 → "bladder" |
| GSE190887 "Urotherlium" cells have wrong condition/technology | Match sci-RNA-seq3 plate barcodes to GEO metadata; set technology = "sci-RNA-seq3" |
| GSE209610 (AKI) condition column | Replace with `condition_level1` |
| GSE119531 (UUO) single sample | Set condition = "UUO_14days" |
| Duplicated samples | Remove: UUO_Day10–14 (also in GSE190887), BladderNormal1/2 |
| Developmental datasets | Remove: E9To13.5Gestation, E18_5_Kidney |
| Spatial lifespan dataset | Remove: MultiOmicSpatialMouseKidney_SexSpecific_lifespan_Chen2025NatGenet |

---

## Sample Categories

Created as `Categories` column in metadata:

| Category | Criteria |
|----------|----------|
| `kidneyMouse` | tissue = kidney, not organoid |
| `bladderMouse` | tissue = bladder, not organoid |
| `ureterMouse` | tissue = ureter, not organoid |
| `bladderOrganoid` | sample_id matches KUDO_GSM3827175/76 |
| `kidneyOrganoid` | sample_id matches KUDO_Yoda/Vehicle/GOF/LOF |
| `ureterOrganoid` | sample_id matches MouseUreterRecon |

---

## Key Parameters

| Parameter | Value | Where set |
|-----------|-------|-----------|
| HVGs (initial integration) | 3 000 | `01`, `05` |
| HVGs (re-integration of gated cells) | 3 000 | `02b` |
| HVGs (scVI export, batch-aware) | 4 000 | `04a` |
| PCA dims (initial) | 20 | `01`, `05` |
| PCA dims (gated re-integration) | 50 | `02b` |
| Harmony batch vars | sample_id + technology (+ source in `01`) | `01`, `02b` |
| Clustering resolution | 0.5 | all |
| scVI latent dims | 30 | `04b` |
| scVI layers | 2 | `04b` |
| scVI gene likelihood | normal (log-norm input) | `04b` |
| DE log2FC threshold | 0.5 | `06`, `07` |
| DE adj. p threshold | 0.05 | `06`, `07` |
| Correlation gene modes | HVG(5000) / allExpr-noMtRibo / allExpr-noMt | `03b` |
| Min expression fraction | 0.10 of samples | `03b` |
| Core gene sharing threshold | ≥ 3 groups | `07` |

---

## Output Files

| File | Produced by | Description |
|------|-------------|-------------|
| `AllUrothelium_harmony_integrated.rds` | `01` | All cells, initial Harmony integration |
| `AllUrothelium_preHarmony_PCA.rds` | `01` | Pre-Harmony checkpoint |
| `AllUrothelium_markers_gated.rds` | `02b` | Urothelium-gated cells, pre-UMAP |
| `AllUrothelium_gated_UMAP.rds` | `05` | Gated cells with UMAP — main downstream input |
| `AllUrothelium_postHarmony_markergated_UMAP.rds` | `02b` | Gated cells UMAP (02b path) |
| `AllUrothelium_gated_UMAP_nosnRNA.rds` | `03b` | scRNA-seq only subset |
| `AllUrothelium_nosnRNA_UMAP.rds` | `02c` | scRNA-only, re-clustered UMAP |
| `AllUrothelium_scvi_input/` | `04a` | MTX + metadata for scVI |
| `AllUrothelium_scvi.h5ad` | `04b` | scVI latent space + UMAP + Leiden clusters |

---

## SLURM Submission Scripts

| Script | Runs |
|--------|------|
| `submit_01_integrate.sh` | `01_integrate_AllUrothelium_harmony.R` |
| `submit_02_markers.sh` | `02b_plot_AllUrothelium_markers_AllscRNA.R` |
| `submit_03_correlations.sh` | `03b_metadata_correlations_AllscRNAUrothelium.R` |
| `submit_04_scvi.sh` | `04a_export_for_scvi.R` → `04b_scvi_AllUrothelium.py` (two-step) |
| `submit_05_HarmonyIntegration.sh` | `05_HarmonyIntegration.R` |

Logs are written to `logs/`.
