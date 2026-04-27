#!/usr/bin/env python3
"""
Extract urothelium-marker-positive cells from scvi_qc_integrated.h5ad.

Markers: Krt5, Krt14, Krt20, Krt8, Krt18, Upk2, Upk3a, Trp63, Upk1a, Upk1b
         (Cd44 is not in the 3000 HVG set and is skipped)

A cell is retained if it has expression > 0 in at least one marker gene.
Output: integration_output/Urothelium_cells.h5ad

Author: Xin Wang
Date:   2026-04-23
"""

import os
import numpy as np
import scanpy as sc
import scipy.sparse as sp

# ── Paths ─────────────────────────────────────────────────────────────────────
base_dir   = "/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
query_path = os.path.join(base_dir, "integration_output", "scvi_qc_integrated.h5ad")
out_path   = os.path.join(base_dir, "integration_output", "Urothelium_cells.h5ad")

# ── Marker genes ──────────────────────────────────────────────────────────────
MARKERS = ["Krt5", "Krt14", "Krt20", "Krt8", "Krt18",
           "Upk2", "Upk3a", "Trp63", "Upk1a", "Upk1b"]
# Cd44 is absent from the 3000 HVG set — skipped

# ── Load ──────────────────────────────────────────────────────────────────────
print("Loading scvi_qc_integrated.h5ad ...")
adata = sc.read_h5ad(query_path)
print(f"  Total cells: {adata.n_obs:,}   genes: {adata.n_vars}")

# ── Check which markers are present ──────────────────────────────────────────
markers_found   = [g for g in MARKERS if g in adata.var_names]
markers_missing = [g for g in MARKERS if g not in adata.var_names]
print(f"  Markers found  ({len(markers_found)}): {markers_found}")
if markers_missing:
    print(f"  Markers absent ({len(markers_missing)}): {markers_missing}")

if not markers_found:
    raise ValueError("None of the marker genes are present in adata.var_names.")

# ── Build positive-cell mask ──────────────────────────────────────────────────
marker_idx = [adata.var_names.get_loc(g) for g in markers_found]
X_markers  = adata.X[:, marker_idx]

if sp.issparse(X_markers):
    positive_mask = np.asarray((X_markers > 0).sum(axis=1)).flatten() > 0
else:
    positive_mask = (X_markers > 0).any(axis=1)

n_positive = positive_mask.sum()
print(f"\n  Cells positive for ≥1 marker: {n_positive:,} / {adata.n_obs:,}"
      f"  ({100 * n_positive / adata.n_obs:.2f}%)")

# ── Per-marker cell counts ────────────────────────────────────────────────────
print("\n  Per-marker positive cell counts:")
for gene, idx in zip(markers_found, marker_idx):
    col = X_markers[:, markers_found.index(gene)]
    if sp.issparse(col):
        n = int((col > 0).sum())
    else:
        n = int((col > 0).sum())
    print(f"    {gene:<10s}: {n:>10,}")

# ── Dataset breakdown ─────────────────────────────────────────────────────────
if "DataSet" in adata.obs.columns:
    ds_counts = adata.obs.loc[positive_mask, "DataSet"].value_counts()
    print("\n  Positive cells per DataSet:")
    print(ds_counts.to_string())

# ── Subset and save ───────────────────────────────────────────────────────────
adata_uro = adata[positive_mask].copy()
print(f"\nSaving {adata_uro.n_obs:,} urothelial cells to:\n  {out_path}")
adata_uro.write_h5ad(out_path)
print("Done.")
