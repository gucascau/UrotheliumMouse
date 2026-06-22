#!/usr/bin/env python3
"""
calc_uro_proportion_h5ad.py

Compute urothelium cell proportions for kidney and bladder datasets.
Reads only the obs metadata from h5ad files via h5py — the expression
matrix is never loaded, so even 27 GB files are processed in seconds.

Outputs (in OUT_DIR):
  uro_proportion_kidney.csv   — per-sample counts and proportions for kidney
  uro_proportion_bladder.csv  — per-sample counts and proportions for bladder
"""

import os
import h5py
import numpy as np
import pandas as pd
from collections import Counter

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR = "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse"

KIDNEY_ALLCELLS = os.path.join(
    BASE_DIR,
    "UsedSingleCells/RenalUrotheliumScripts/output/RenalUrothelium_allcells_scvi.h5ad"
)
KIDNEY_URO = os.path.join(
    BASE_DIR,
    "UsedSingleCells/RenalUrotheliumScripts/output/RenalUrothelium_uro_cells_fullgene_scvi.h5ad"
)
BLADDER_ALLCELLS = os.path.join(
    BASE_DIR,
    "UsedSingleCellsRawResults/BladderUrothelium/BladderUrotheliumScripts/output/"
    "BladderUrothelium_allcells_scvi.h5ad"
)
BLADDER_URO = os.path.join(
    BASE_DIR,
    "FinalUrotheliumCells/BladderUrothelium_uro_cells_scvi.h5ad"
)

OUT_DIR = os.path.join(
    BASE_DIR,
    "UsedSingleCells/UrotheliumIntegrationScripts/output/uro_proportion"
)
os.makedirs(OUT_DIR, exist_ok=True)


# ── h5ad obs reader ───────────────────────────────────────────────────────────

def read_obs_col(h5file, col):
    """Return a pandas array for one obs column (categorical or string)."""
    grp = h5file["obs"][col]
    if "codes" in grp:
        cats = list(grp["categories"].asstr())
        codes = np.array(grp["codes"])
        return pd.Categorical.from_codes(codes, cats)
    return list(grp.asstr())


def read_obs(h5_path, cols):
    """
    Read selected obs columns from an h5ad without loading X.
    Returns a DataFrame indexed by cell barcodes.
    """
    with h5py.File(h5_path, "r") as f:
        obs_keys = set(f["obs"].keys())

        # Cell barcodes (index)
        idx_key = "_index" if "_index" in obs_keys else sorted(obs_keys)[0]
        index = read_obs_col(f, idx_key)

        data = {}
        for col in cols:
            if col in obs_keys:
                data[col] = read_obs_col(f, col)
            else:
                data[col] = None

    return pd.DataFrame(data, index=index)


# ── Core computation ──────────────────────────────────────────────────────────

META_COLS = ["orig.ident", "gsm_id", "condition", "technology", "tissue", "paper", "source_GEO"]

def compute_proportion(label, allcells_path, uro_path, sample_col="sample_id"):
    """
    For each sample_id in allcells_path, report:
      total_cells, uro_cells, uro_proportion.
    Includes extra metadata columns when present.
    """
    print(f"\n{'='*60}")
    print(f"  {label}")
    print(f"  allcells : {os.path.basename(allcells_path)}")
    print(f"  uro      : {os.path.basename(uro_path)}")
    print(f"{'='*60}")

    for p in (allcells_path, uro_path):
        if not os.path.exists(p):
            raise FileNotFoundError(f"File not found: {p}")

    print("  Reading allcells obs ...")
    all_obs = read_obs(allcells_path, [sample_col] + META_COLS)
    print(f"  → {len(all_obs):,} total cells across "
          f"{all_obs[sample_col].nunique()} samples")

    print("  Reading urothelium obs ...")
    uro_obs = read_obs(uro_path, [sample_col] + META_COLS)
    print(f"  → {len(uro_obs):,} urothelium cells across "
          f"{uro_obs[sample_col].nunique()} samples")

    # Cell counts
    total_cnt = all_obs.groupby(sample_col, observed=True).size().rename("total_cells")
    uro_cnt   = uro_obs.groupby(sample_col, observed=True).size().rename("uro_cells")

    # Best available metadata per sample (first non-null value)
    keep_meta = [c for c in META_COLS if c in all_obs.columns and all_obs[c].notna().any()]
    meta = (
        all_obs.groupby(sample_col, observed=True)[keep_meta]
        .first()
        .reset_index()
        .set_index(sample_col)
    )

    df = pd.concat([total_cnt, uro_cnt, meta], axis=1)
    df["uro_cells"]      = df["uro_cells"].fillna(0).astype(int)
    df["uro_proportion"] = (df["uro_cells"] / df["total_cells"]).round(6)
    df["uro_pct"]        = (df["uro_proportion"] * 100).round(4)

    df.insert(0, "dataset_type", label)
    df.index.name = "sample_id"
    df = df.reset_index().sort_values("sample_id")

    # Print summary
    print(f"\n  {'sample_id':<50} {'total':>8} {'uro':>8} {'pct':>8}")
    print(f"  {'-'*50} {'-'*8} {'-'*8} {'-'*8}")
    for _, row in df.iterrows():
        print(f"  {row['sample_id']:<50} {row['total_cells']:>8,} "
              f"{row['uro_cells']:>8,} {row['uro_pct']:>7.2f}%")
    print(f"  {'TOTAL':<50} {df['total_cells'].sum():>8,} "
          f"{df['uro_cells'].sum():>8,} "
          f"{100*df['uro_cells'].sum()/df['total_cells'].sum():>7.2f}%")

    return df


# ── Run ───────────────────────────────────────────────────────────────────────

print("=== Urothelium Proportion Calculator (h5ad) ===")

kidney_df = compute_proportion("kidney", KIDNEY_ALLCELLS, KIDNEY_URO)
kidney_out = os.path.join(OUT_DIR, "uro_proportion_kidney.csv")
kidney_df.to_csv(kidney_out, index=False)
print(f"\n  Saved: {kidney_out}")

bladder_df = compute_proportion("bladder", BLADDER_ALLCELLS, BLADDER_URO)
bladder_out = os.path.join(OUT_DIR, "uro_proportion_bladder.csv")
bladder_df.to_csv(bladder_out, index=False)
print(f"\n  Saved: {bladder_out}")

print("\n=== Done ===")
