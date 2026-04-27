#!/usr/bin/env python3
"""
03_extract_uro.py — Filter the scVI-integrated AnnData to cells that
express at least one urothelial marker gene (expression > 0).

Markers: Krt5, Krt14, Krt20, Krt8, Krt18, Upk2, Upk3a, Trp63, Upk1a, Upk1b

Input  : output/RenalUrothelium_integrated.h5ad
Output : output/RenalUrothelium_uro_cells.h5ad
"""

import os
import numpy as np
import scanpy as sc
import scipy.sparse as sp

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR   = "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
SCRIPT_DIR = os.path.join(BASE_DIR, "RenalUrotheliumScripts")
OUT_DIR    = os.path.join(SCRIPT_DIR, "output")

IN_PATH  = os.path.join(OUT_DIR, "RenalUrothelium_integrated.h5ad")
OUT_PATH = os.path.join(OUT_DIR, "RenalUrothelium_uro_cells.h5ad")

# ── Markers ───────────────────────────────────────────────────────────────────
MARKERS = ["Krt5", "Krt14", "Krt20", "Krt8", "Krt18",
           "Upk2", "Upk3a", "Trp63", "Upk1a", "Upk1b"]

###############################################################################
# Load
###############################################################################

print(f"Loading {IN_PATH} ...")
adata = sc.read_h5ad(IN_PATH)
print(f"  Total: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

###############################################################################
# Marker check
###############################################################################

found   = [g for g in MARKERS if g in adata.var_names]
missing = [g for g in MARKERS if g not in adata.var_names]
print(f"\n  Markers found  ({len(found)}): {found}")
if missing:
    print(f"  Markers absent ({len(missing)}): {missing}")

if not found:
    raise ValueError("None of the urothelium markers are in adata.var_names.")

###############################################################################
# Build positive-cell mask
###############################################################################

idx    = [adata.var_names.get_loc(g) for g in found]
X_mark = adata.X[:, idx]

if sp.issparse(X_mark):
    mask = np.asarray((X_mark > 0).sum(axis=1)).flatten() > 0
else:
    mask = (X_mark > 0).any(axis=1)

n = mask.sum()
print(f"\n  Urothelium-positive: {n:,} / {adata.n_obs:,} "
      f"({100 * n / adata.n_obs:.2f}%)")

print("\n  Per-marker positive counts:")
for gene, i in zip(found, idx):
    col = X_mark[:, found.index(gene)]
    nc  = int((col > 0).sum()) if not sp.issparse(col) else int((col > 0).nnz if hasattr(col,'nnz') else (col > 0).sum())
    print(f"    {gene:<10s}: {nc:>8,}")

print("\n  Positive cells per sample:")
print(adata.obs.loc[mask, "sample_id"].value_counts().to_string())

if "condition" in adata.obs.columns:
    print("\n  Positive cells per condition:")
    print(adata.obs.loc[mask, "condition"].value_counts().to_string())

if "leiden_scVI" in adata.obs.columns:
    print("\n  Positive cells per scVI cluster:")
    print(adata.obs.loc[mask, "leiden_scVI"].value_counts().sort_index().to_string())

###############################################################################
# Subset and save
###############################################################################

adata_uro = adata[mask].copy()
print(f"\nSaving {adata_uro.n_obs:,} urothelial cells → {OUT_PATH}")
adata_uro.write_h5ad(OUT_PATH)
print("Done.")
