#!/usr/bin/env python3
"""
01_preprocess.py — Load all h5ad files from RenalUrothelium/, harmonise
metadata, extract raw integer counts, select batch-aware HVGs, and save
a single concatenated h5ad ready for scVI training.

All 19 datasets are included.  Sources and their count strategies:

  KidneyHealthy1-5, KidneyTET2UUO, KidneyUUO1-8, KidneyrUUO1
      → integer counts in X  (standard 10X h5ad)

  MouseKidneyATLAS_MKA_updated.h5ad      (141k cells, 16k genes)
  mouse_kidney_snRNAseq_Lake2025_*.h5ad  (316k cells, 28k genes)
  MultiOmicSpatialMouseKidney_*.h5ad     (203k cells, 31k genes)
      → CellXGene format: X = normalised; raw.X = integer counts (float32)
        var_names = Ensembl IDs → remapped to gene symbols via feature_name

  UUOProjectObject_F.h5ad                (310k cells, 18k genes)
      → Seurat export: X = integer counts (float64), sample from orig.ident

Outputs
-------
  output/RenalUrothelium_preprocessed.h5ad
"""

import gc
import os
import re
import glob
import numpy as np
import anndata as ad
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR = "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
DATA_DIR = os.path.join(BASE_DIR, "RenalUrothelium")
OUT_DIR  = os.path.join(BASE_DIR, "RenalUrotheliumScripts", "output")
os.makedirs(OUT_DIR, exist_ok=True)

OUT_PATH = os.path.join(OUT_DIR, "RenalUrothelium_preprocessed.h5ad")

# ── Parameters ────────────────────────────────────────────────────────────────
N_HVG     = 3000
BATCH_KEY = "sample_id"

# ── Per-filename metadata (condition, technology) ─────────────────────────────
SAMPLE_META = [
    (r"KidneyHealthy\d+",   "Healthy",  "10X"),
    (r"KidneyTET2UUO",      "TET2UUO",  "10X"),
    (r"KidneyUUO\d+",       "UUO",      "10X"),
    (r"KidneyrUUO\d+",      "rUUO",     "10X"),
    (r".*MKA.*",            "Healthy",  "10X"),
    (r".*Lake.*",           "Healthy",  "snRNA-seq"),
    (r".*Chen.*",           "Healthy",  "Spatial"),
    (r"UUOProjectObject.*", "UUO",      "10X"),
]

def infer_meta(sid):
    for pattern, condition, technology in SAMPLE_META:
        if re.fullmatch(pattern, sid, re.IGNORECASE):
            return condition, technology
    return "Unknown", "10X"


# ── Helper: remap Ensembl IDs → gene symbols (CellXGene format) ──────────────
def use_gene_symbols(adata, sid):
    """If var has a 'feature_name' column, set it as var_names (gene symbols)."""
    if "feature_name" in adata.var.columns:
        adata.var_names = adata.var["feature_name"].astype(str).values
        adata.var_names_make_unique()
        print(f"  [{sid}] var_names remapped to gene symbols via feature_name")
    return adata


# ── Helper: extract raw integer counts ───────────────────────────────────────
def prepare_counts(adata, sid):
    """
    Populate adata.layers['counts'] with raw integer counts and set X to
    log-normalised values for HVG selection.

    Priority:
      1. adata.raw.X   — CellXGene datasets store raw counts here
      2. existing 'counts' layer with integer values
      3. X is already integer (standard 10X export)
      4. X is float (pre-normalised) — back-transform to pseudo-counts
    """
    # 1. Use adata.raw if available (CellXGene format: raw has full gene × int counts)
    if adata.raw is not None:
        raw_adata = adata.raw.to_adata()
        use_gene_symbols(raw_adata, sid + "/raw")

        raw_x = raw_adata.X
        vals  = raw_x.data if sp.issparse(raw_x) else raw_x.ravel()
        vals  = vals[vals != 0][:5000]
        if len(vals) > 0 and np.all(vals == np.floor(vals)):
            print(f"  [{sid}] using raw.X (integer counts, "
                  f"{raw_adata.n_vars:,} genes)")
            # Rebuild adata from raw to get full gene set with log-norm X
            counts = sp.csr_matrix(raw_x)
            adata_new = ad.AnnData(
                X   = counts.copy().astype(np.float32),
                obs = adata.obs.copy(),
                var = raw_adata.var.copy(),
            )
            adata_new.layers["counts"] = counts
            # Normalise for HVG / scVI log-norm layer
            sc.pp.normalize_total(adata_new, target_sum=1e4)
            sc.pp.log1p(adata_new)
            return adata_new

    # 2. Existing integer counts layer
    if "counts" in adata.layers:
        lyr  = adata.layers["counts"]
        vals = lyr.data if sp.issparse(lyr) else lyr.ravel()
        vals = vals[vals != 0][:5000]
        if len(vals) > 0 and np.all(vals == np.floor(vals)):
            print(f"  [{sid}] using existing integer 'counts' layer")
            adata.X = adata.layers["counts"].copy().astype(np.float32)
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            return adata

    # 3. X is raw integer counts
    x    = adata.X
    vals = x.data if sp.issparse(x) else x.ravel()
    vals = vals[vals != 0][:5000]
    if len(vals) > 0 and np.all(vals == np.floor(vals)):
        print(f"  [{sid}] raw integer counts in X")
        counts = (sp.csr_matrix(x) if not sp.issparse(x) else x).astype(np.float32)
        adata.X = None  # free original before copying
        adata.layers["counts"] = counts
        adata.X = counts.copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        return adata

    # 4. Float X — create pseudo-counts
    x_max = float(vals.max()) if len(vals) > 0 else 0
    if x_max < 20:
        print(f"  [{sid}] log-norm X (max={x_max:.2f}) → pseudo-counts via expm1×10")
        x_arr   = x.toarray() if sp.issparse(x) else np.array(x, dtype=np.float32)
        pseudo  = np.round(np.clip(np.expm1(x_arr) * 10, 0, None)).astype(np.float32)
    else:
        print(f"  [{sid}] large-scale normalised X (max={x_max:.2f}) → round as pseudo-counts")
        x_arr  = x.toarray() if sp.issparse(x) else np.array(x, dtype=np.float32)
        pseudo = np.round(np.clip(x_arr, 0, None)).astype(np.float32)

    adata.layers["counts"] = sp.csr_matrix(pseudo)
    print(f"  [{sid}] X retained as log-norm for HVG selection")
    return adata


# ── Helper: harmonise obs metadata ───────────────────────────────────────────
def harmonise_obs(adata, sid, condition, technology):
    """Assign standardised sample_id / condition / technology / cell_type_original."""

    # sample_id: for UUOProjectObject use orig.ident (per-sample identity)
    if "orig.ident" in adata.obs.columns:
        adata.obs["sample_id"] = adata.obs["orig.ident"].astype(str)
    else:
        adata.obs["sample_id"] = sid

    # technology: prefer existing obs column if informative
    if "technology" in adata.obs.columns:
        adata.obs["technology"] = adata.obs["technology"].astype(str)
    else:
        adata.obs["technology"] = technology

    # condition: prefer existing obs column if informative
    if "condition" in adata.obs.columns:
        adata.obs["condition"] = adata.obs["condition"].astype(str)
    else:
        adata.obs["condition"] = condition

    # source: the filename stem (useful as a batch covariate for reference datasets)
    adata.obs["source"] = sid

    # pct_mt: carry over if present (Seurat export uses percent.mt)
    if "percent.mt" in adata.obs.columns and "pct_mt" not in adata.obs.columns:
        adata.obs["pct_mt"] = adata.obs["percent.mt"]

    # cell_type_original: unify annotation columns across all datasets
    if "cell_type_original" not in adata.obs.columns:
        for col in ["author_cell_type", "celltype_final", "CellStateLevel1",
                    "cell_type", "CellType", "ident", "leiden"]:
            if col in adata.obs.columns:
                adata.obs["cell_type_original"] = adata.obs[col].astype(str)
                print(f"  [{sid}] cell_type_original ← {col}")
                break

    return adata


def _coerce_is_primary_data(series):
    """Convert CellxGene-style mixed/object values into nullable booleans."""
    mapping = {
        True: True,
        False: False,
        1: True,
        0: False,
        "1": True,
        "0": False,
        "true": True,
        "false": False,
        "True": True,
        "False": False,
        "": pd.NA,
        "nan": pd.NA,
        "None": pd.NA,
        "none": pd.NA,
    }
    coerced = series.map(lambda x: mapping.get(x, pd.NA))
    return coerced.astype("boolean")


def sanitise_dataframe_for_h5ad(df, axis_name):
    """Rename reserved columns and coerce object metadata into HDF5-safe dtypes."""
    df = df.copy()

    if "_index" in df.columns:
        replacement = f"{axis_name}_index_original"
        suffix = 1
        while replacement in df.columns:
            replacement = f"{axis_name}_index_original_{suffix}"
            suffix += 1
        df = df.rename(columns={"_index": replacement})
        print(f"  [{axis_name}] renamed reserved column '_index' → '{replacement}'")

    for col in df.columns:
        series = df[col]
        if col == "is_primary_data":
            df[col] = _coerce_is_primary_data(series)
            continue
        if series.dtype == "object":
            non_na = series.dropna()
            if non_na.empty:
                df[col] = series.fillna("").astype(str)
                continue

            if non_na.map(lambda x: isinstance(x, (bool, np.bool_))).all():
                df[col] = series.astype("boolean")
                continue

            if non_na.map(lambda x: isinstance(x, (int, float, np.integer, np.floating, bool, np.bool_))).all():
                df[col] = pd.to_numeric(series, errors="coerce")
                continue

            df[col] = series.fillna("").map(str)

    return df


def sanitise_anndata_for_h5ad(adata):
    """Normalise obs/var metadata so AnnData can be written reliably."""
    adata.obs = sanitise_dataframe_for_h5ad(adata.obs, "obs")
    adata.var = sanitise_dataframe_for_h5ad(adata.var, "var")
    return adata


###############################################################################
# STEP 1: Load and preprocess each file
###############################################################################

h5ad_files = sorted(glob.glob(os.path.join(DATA_DIR, "*.h5ad")))
print(f"Found {len(h5ad_files)} h5ad files in {DATA_DIR}\n")

adatas  = []

for path in h5ad_files:
    fname = os.path.basename(path)
    sid   = os.path.splitext(fname)[0]

    print(f"{'='*60}")
    print(f"Loading {fname} ...")
    a = sc.read_h5ad(path)
    print(f"  {a.n_obs:,} cells × {a.n_vars:,} genes  | "
          f"layers={list(a.layers.keys())} | "
          f"raw={'yes' if a.raw is not None else 'no'}")

    # Remap var_names to gene symbols for CellXGene-format files
    a = use_gene_symbols(a, sid)

    # Prepare counts layer + log-norm X
    a = prepare_counts(a, sid)

    # Harmonise metadata
    condition, technology = infer_meta(sid)
    a = harmonise_obs(a, sid, condition, technology)

    print(f"  [{sid}] ready | {a.n_obs:,} cells × {a.n_vars:,} genes | "
          f"condition={condition} | technology={technology}")
    adatas.append(a)
    gc.collect()

print(f"\nTotal samples loaded: {len(adatas)}")


###############################################################################
# STEP 2: Concatenate (outer join — missing genes in any dataset filled with 0)
###############################################################################

print("\nConcatenating ...")
adata = ad.concat(
    adatas,
    join         = "outer",
    merge        = "first",
    index_unique = "-",      # makes cell barcodes unique across batches
)
del adatas
gc.collect()

# Fill NaN introduced by outer join with 0
for layer in list(adata.layers.keys()):
    lyr = adata.layers[layer]
    if sp.issparse(lyr):
        lyr.data = np.nan_to_num(lyr.data, nan=0.0)
    else:
        adata.layers[layer] = np.nan_to_num(lyr, nan=0.0)
if sp.issparse(adata.X):
    adata.X.data = np.nan_to_num(adata.X.data, nan=0.0)
else:
    adata.X = np.nan_to_num(adata.X, nan=0.0)

print(f"  Concatenated: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
print("\n  Cells per source:")
print(adata.obs["source"].value_counts().to_string())
print("\n  Cells per condition:")
print(adata.obs["condition"].value_counts().to_string())
print("\n  Cells per technology:")
print(adata.obs["technology"].value_counts().to_string())


###############################################################################
# STEP 3: Batch-aware HVG selection
###############################################################################

print(f"\nSelecting {N_HVG} HVGs (batch_key={BATCH_KEY}, flavour=seurat_v3) ...")
sc.pp.highly_variable_genes(
    adata,
    n_top_genes = N_HVG,
    batch_key   = BATCH_KEY,
    flavor      = "seurat_v3",
    layer       = "counts",
    subset      = False,
)
n_hvg = int(adata.var["highly_variable"].sum())
print(f"  HVGs selected: {n_hvg}")

adata = adata[:, adata.var["highly_variable"]].copy()
print(f"  After HVG subset: {adata.n_obs:,} cells × {adata.n_vars:,} genes")


###############################################################################
# STEP 4: Save
###############################################################################

adata = sanitise_anndata_for_h5ad(adata)
print(f"\nSaving to {OUT_PATH} ...")
adata.write_h5ad(OUT_PATH)
print("Done.")
print(f"\n===== Preprocessing complete =====")
print(f"  Cells : {adata.n_obs:,}")
print(f"  Genes : {adata.n_vars:,}")
print(f"  Output: {OUT_PATH}")
