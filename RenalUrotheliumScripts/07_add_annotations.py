#!/usr/bin/env python3
"""
07_add_annotations.py — Transfer scANVI annotations from Lake and MKA
reference models into RenalUrothelium_integrated.h5ad.

Inputs:
  output/RenalUrothelium_integrated.h5ad
  output/annotation_Lake/RenalUrothelium_scanvi_Lake_annotations.csv
  output/annotation_MKA/RenalUrothelium_scanvi_MKA_annotations.csv

New obs columns added:
  scanvi_Lake_label      : predicted cell type from Lake reference
  scanvi_Lake_confidence : prediction confidence (0–1)
  scanvi_MKA_label       : predicted cell type from MKA reference
  scanvi_MKA_confidence  : prediction confidence (0–1)

Output:
  output/RenalUrothelium_integrated.h5ad  (overwritten in-place)
"""

import os
import pandas as pd
import scanpy as sc

BASE_DIR   = "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
SCRIPT_DIR = os.path.join(BASE_DIR, "RenalUrotheliumScripts")
OUT_DIR    = os.path.join(SCRIPT_DIR, "output")

H5AD_IN    = os.path.join(OUT_DIR, "RenalUrothelium_integrated.h5ad")
H5AD_OUT   = os.path.join(OUT_DIR, "RenalUrothelium_integrated_annotated.h5ad")
LAKE_CSV   = os.path.join(OUT_DIR, "annotation_Lake", "RenalUrothelium_scanvi_Lake_annotations.csv")
MKA_CSV    = os.path.join(OUT_DIR, "annotation_MKA",  "RenalUrothelium_scanvi_MKA_annotations.csv")

# ── Load integrated AnnData ───────────────────────────────────────────────────
print(f"Loading {H5AD_IN} ...")
adata = sc.read_h5ad(H5AD_IN)
print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")
print(f"  Existing obs columns: {list(adata.obs.columns)}")

# ── Load annotation CSVs (index = cell barcode) ───────────────────────────────
lake = pd.read_csv(LAKE_CSV, index_col=0)[["scanvi_Lake_label", "scanvi_Lake_confidence"]]
mka  = pd.read_csv(MKA_CSV,  index_col=0)[["scanvi_MKA_label",  "scanvi_MKA_confidence"]]

# The annotation scripts append '-query' to barcodes; strip it to match adata.obs_names
lake.index = lake.index.str.replace(r"-query$", "", regex=True)
mka.index  = mka.index.str.replace(r"-query$", "", regex=True)

print(f"\nLake annotations: {len(lake):,} rows")
print(f"  Label counts:\n{lake['scanvi_Lake_label'].value_counts().head(10).to_string()}")

print(f"\nMKA annotations: {len(mka):,} rows")
print(f"  Label counts:\n{mka['scanvi_MKA_label'].value_counts().head(10).to_string()}")

# ── Join onto adata.obs (left join — cells not in CSV get NaN) ────────────────
for col in ["scanvi_Lake_label", "scanvi_Lake_confidence",
            "scanvi_MKA_label",  "scanvi_MKA_confidence"]:
    if col in adata.obs.columns:
        adata.obs.drop(columns=col, inplace=True)

adata.obs = adata.obs.join(lake, how="left")
adata.obs = adata.obs.join(mka,  how="left")

# Report coverage
n_lake = adata.obs["scanvi_Lake_label"].notna().sum()
n_mka  = adata.obs["scanvi_MKA_label"].notna().sum()
print(f"\nCoverage after join:")
print(f"  Lake annotation: {n_lake:,} / {adata.n_obs:,} cells "
      f"({100 * n_lake / adata.n_obs:.1f}%)")
print(f"  MKA  annotation: {n_mka:,} / {adata.n_obs:,} cells "
      f"({100 * n_mka / adata.n_obs:.1f}%)")

print(f"\nTop Lake labels across all cells:")
print(adata.obs["scanvi_Lake_label"].value_counts().head(15).to_string())

print(f"\nTop MKA labels across all cells:")
print(adata.obs["scanvi_MKA_label"].value_counts().head(15).to_string())

# ── Save ──────────────────────────────────────────────────────────────────────
print(f"\nSaving annotated AnnData → {H5AD_OUT} ...")
adata.write_h5ad(H5AD_OUT)
print("Done.")
