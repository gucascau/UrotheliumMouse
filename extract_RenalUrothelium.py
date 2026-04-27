#!/usr/bin/env python3
"""
Extract renal urothelial cells for Harmony integration.

Sources
-------
Kidney scVI cells (from scvi_qc_integrated.h5ad, X = log-norm, 3000 HVGs):
  KidneyHealthy1-5, KidneyTET2UUO, KidneyUUO1-6, KidneyrUUO1

Reference datasets (from qc_h5ad/, pre-normalised, all genes):
  LakesnRNA, ChenSpatial, MKA

Workflow per source
-------------------
1. Load AnnData
2. Filter to urothelium-marker-positive cells (≥1 marker > 0)
3. Normalise if raw integer counts (normalize_total + log1p); otherwise use X as-is
4. Subset genes to the 3000-gene HVG list from the scVI h5ad
5. Ensure consistent metadata columns (sample_id, technology, condition, source)

Then concatenate all sources (outer join; missing genes → 0) and save.

Output: integration_output/RenalUrothelium_cells.h5ad

Author: Xin Wang
Date:   2026-04-23
"""

import os
import numpy as np
import anndata as ad
import scanpy as sc
import scipy.sparse as sp

# ── Paths ─────────────────────────────────────────────────────────────────────
base_dir   = "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
scvi_path  = os.path.join(base_dir, "integration_output", "scvi_qc_integrated.h5ad")
qc_dir     = os.path.join(base_dir, "qc_h5ad")
out_path   = os.path.join(base_dir, "integration_output", "RenalUrothelium_cells.h5ad")

# ── Target kidney samples ─────────────────────────────────────────────────────
TARGET_SAMPLES = [
    "KidneyHealthy1", "KidneyHealthy2", "KidneyHealthy3",
    "KidneyHealthy4", "KidneyHealthy5",
    "KidneyTET2UUO",
    "KidneyUUO1", "KidneyUUO2", "KidneyUUO3",
    "KidneyUUO4", "KidneyUUO5", "KidneyUUO6",
    "KidneyrUUO1",
]

# ── Reference datasets (sample_id, technology, condition) ─────────────────────
REF_DATASETS = {
    "LakesnRNA":  {"path": os.path.join(qc_dir, "LakesnRNA.h5ad"),
                   "technology": "snRNA-seq", "condition": "Healthy"},
    "ChenSpatial":{"path": os.path.join(qc_dir, "ChenSpatial.h5ad"),
                   "technology": "Spatial",   "condition": "Healthy"},
    "MKA":        {"path": os.path.join(qc_dir, "MKA.h5ad"),
                   "technology": "10X",        "condition": "Healthy"},
}

# ── Urothelium markers ────────────────────────────────────────────────────────
MARKERS = ["Krt5", "Krt14", "Krt20", "Krt8", "Krt18",
           "Upk2", "Upk3a", "Trp63", "Upk1a", "Upk1b"]


# ── Helper: filter to urothelium-positive cells ───────────────────────────────
def filter_uro(adata, label):
    found   = [g for g in MARKERS if g in adata.var_names]
    missing = [g for g in MARKERS if g not in adata.var_names]
    print(f"  [{label}] markers found ({len(found)}): {found}")
    if missing:
        print(f"  [{label}] markers absent: {missing}")
    if not found:
        print(f"  [{label}] WARNING — no markers present, skipping.")
        return None

    idx    = [adata.var_names.get_loc(g) for g in found]
    X_mark = adata.X[:, idx]
    if sp.issparse(X_mark):
        mask = np.asarray((X_mark > 0).sum(axis=1)).flatten() > 0
    else:
        mask = (X_mark > 0).any(axis=1)

    n = mask.sum()
    print(f"  [{label}] urothelium-positive: {n:,} / {adata.n_obs:,} "
          f"({100 * n / adata.n_obs:.2f}%)")
    return adata[mask].copy()


# ── Helper: normalise if raw counts ──────────────────────────────────────────
def ensure_lognorm(adata, label):
    x_data = adata.X.data if sp.issparse(adata.X) else adata.X.ravel()
    x_data = x_data[x_data != 0]
    if len(x_data) > 0 and np.all(x_data == np.floor(x_data)):
        print(f"  [{label}] raw counts detected — running normalize_total + log1p")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    else:
        print(f"  [{label}] pre-normalised data — using X as-is")
    return adata


# ── Helper: enforce required metadata columns ─────────────────────────────────
REQUIRED_COLS = ["sample_id", "technology", "condition", "source"]

def harmonise_obs(adata, sample_id, technology, condition, source):
    adata.obs["sample_id"]  = sample_id
    adata.obs["technology"] = technology
    adata.obs["condition"]  = condition
    adata.obs["source"]     = source
    # Keep only a clean subset of obs columns to avoid concat conflicts
    keep = [c for c in REQUIRED_COLS if c in adata.obs.columns]
    # Also preserve annotation columns if present
    for extra in ["cell_type_original", "cell_type", "leiden", "CellType"]:
        if extra in adata.obs.columns:
            keep.append(extra)
    adata.obs = adata.obs[keep].copy()
    return adata


###############################################################################
# PART 1: Kidney cells from scvi_qc_integrated.h5ad
###############################################################################

print("=" * 70)
print("PART 1 — Kidney scVI cells")
print("=" * 70)

print(f"Loading {scvi_path} ...")
adata_scvi = sc.read_h5ad(scvi_path)
print(f"  Total cells: {adata_scvi.n_obs:,}   genes: {adata_scvi.n_vars}")

# Record the 3000-gene HVG list — all reference datasets will be subset to this
HVG_GENES = list(adata_scvi.var_names)
print(f"  HVG gene list: {len(HVG_GENES)} genes")

# Filter to target kidney samples
found_samples   = [s for s in TARGET_SAMPLES if s in adata_scvi.obs["sample_id"].values]
missing_samples = [s for s in TARGET_SAMPLES if s not in adata_scvi.obs["sample_id"].values]
print(f"\n  Samples present ({len(found_samples)}/{len(TARGET_SAMPLES)}): {found_samples}")
if missing_samples:
    print(f"  WARNING — not found: {missing_samples}")

adata_scvi = adata_scvi[adata_scvi.obs["sample_id"].isin(found_samples)].copy()
print(f"  After sample filter: {adata_scvi.n_obs:,} cells")

# Filter to urothelium-positive cells
adata_kidney = filter_uro(adata_scvi, "scVI-kidney")
if adata_kidney is None:
    raise ValueError("No urothelium cells found in kidney scVI data.")
del adata_scvi

print("\n  Cells per sample:")
print(adata_kidney.obs["sample_id"].value_counts().to_string())

# Harmonise metadata — carry condition/technology from existing obs if present
for col in ["technology", "condition"]:
    if col not in adata_kidney.obs.columns:
        adata_kidney.obs[col] = "Unknown"

adata_kidney = harmonise_obs(
    adata_kidney,
    sample_id  = adata_kidney.obs["sample_id"],   # keep per-cell value
    technology = adata_kidney.obs["technology"],
    condition  = adata_kidney.obs["condition"],
    source     = "scVI",
)
print(f"\n  Kidney uro cells retained: {adata_kidney.n_obs:,}")


###############################################################################
# PART 2: Reference datasets (LakesnRNA, ChenSpatial, MKA)
###############################################################################

print("\n" + "=" * 70)
print("PART 2 — Reference datasets")
print("=" * 70)

ref_adatas = []

for sid, info in REF_DATASETS.items():
    path = info["path"]
    print(f"\n--- {sid} ---")
    if not os.path.exists(path):
        print(f"  WARNING — file not found, skipping: {path}")
        continue

    ref = sc.read_h5ad(path)
    print(f"  Loaded: {ref.n_obs:,} cells × {ref.n_vars} genes")

    # Normalise if needed
    ref = ensure_lognorm(ref, sid)

    # Filter to urothelium-positive cells
    ref_uro = filter_uro(ref, sid)
    del ref
    if ref_uro is None:
        continue

    # Subset genes to the 3000 HVG gene list
    common_genes = [g for g in HVG_GENES if g in ref_uro.var_names]
    print(f"  [{sid}] genes overlapping HVG list: {len(common_genes):,} / {len(HVG_GENES)}")
    if len(common_genes) == 0:
        print(f"  [{sid}] WARNING — no gene overlap, skipping.")
        continue
    ref_uro = ref_uro[:, common_genes].copy()

    # Harmonise metadata
    ref_uro = harmonise_obs(
        ref_uro,
        sample_id  = sid,
        technology = info["technology"],
        condition  = info["condition"],
        source     = sid,
    )
    print(f"  [{sid}] urothelium cells kept: {ref_uro.n_obs:,}")
    ref_adatas.append(ref_uro)


###############################################################################
# PART 3: Concatenate and save
###############################################################################

print("\n" + "=" * 70)
print("PART 3 — Concatenate & save")
print("=" * 70)

all_adatas = [adata_kidney] + ref_adatas

print(f"\nConcatenating {len(all_adatas)} batches (outer join on genes) ...")
merged = ad.concat(
    all_adatas,
    join         = "outer",   # missing genes in references → 0
    merge        = "first",
    index_unique = "-",
)

# Fill NaNs introduced by outer join
if sp.issparse(merged.X):
    merged.X.data = np.nan_to_num(merged.X.data)
else:
    merged.X = np.nan_to_num(merged.X)

print(f"  Merged: {merged.n_obs:,} cells × {merged.n_vars} genes")
print("\n  Cells per source:")
print(merged.obs["source"].value_counts().to_string())
print("\n  Cells per sample_id:")
print(merged.obs["sample_id"].value_counts().to_string())
if "condition" in merged.obs.columns:
    print("\n  Cells per condition:")
    print(merged.obs["condition"].value_counts().to_string())

print(f"\nSaving to:\n  {out_path}")
merged.write_h5ad(out_path)
print("Done.")
