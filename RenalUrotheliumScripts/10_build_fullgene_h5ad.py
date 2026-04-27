#!/usr/bin/env python3
"""
10_build_fullgene_h5ad.py — Build a full-transcriptome (~25k gene) h5ad for
experimental cells only, without re-running the integration pipeline.

Strategy
--------
  1. Load obs + reductions from the integrated annotated h5ad (backed, no X).
  2. Keep only the 14 experimental samples (exclude reference atlases).
  3. For each sample, fetch raw counts for ALL genes from the original h5ad.
  4. Concatenate with outer join (missing genes in any sample → 0).
  5. Normalise per cell (normalize_total + log1p) → store in layers["lognorm"].
  6. Transfer metadata and UMAP embeddings from the integrated object.
  7. Save as output/RenalUrothelium_experimental_fullgene.h5ad.

This h5ad can then be converted to a Seurat RDS by 10b_convert_fullgene_to_rds.R
which also transfers the scVI/scANVI UMAP reductions.

Output
------
  output/RenalUrothelium_experimental_fullgene.h5ad
    layers["counts"]   — raw integer counts (~25k genes)
    layers["lognorm"]  — log-normalised (normalize_total 10k + log1p)
    X                  — same as lognorm (scanpy default)
    obs                — all metadata from integrated object
    obsm               — UMAP embeddings from integrated object
"""

import gc
import os
import sys
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scipy.sparse as sp

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR   = "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
DATA_DIR   = os.path.join(BASE_DIR, "RenalUrothelium")
SCRIPT_DIR = os.path.join(BASE_DIR, "RenalUrotheliumScripts")
OUT_DIR    = os.path.join(SCRIPT_DIR, "output")
ANNOT_H5AD = os.path.join(OUT_DIR, "RenalUrothelium_integrated_annotated.h5ad")
OUT_PATH   = os.path.join(OUT_DIR, "RenalUrothelium_experimental_fullgene.h5ad")

# ── Experimental samples ───────────────────────────────────────────────────────
EXPERIMENTAL_SOURCES = {
    "KidneyHealthy1", "KidneyHealthy2", "KidneyHealthy3",
    "KidneyHealthy4", "KidneyHealthy5",
    "KidneyTET2UUO",
    "KidneyUUO1",  "KidneyUUO2",  "KidneyUUO3",  "KidneyUUO4",
    "KidneyUUO5",  "KidneyUUO6",  "KidneyUUO7",  "KidneyUUO8",
    "KidneyrUUO1",
}

# UMAP and reduction keys to carry over from the integrated object
OBSM_KEYS = ["X_scVI", "X_scANVI", "X_umap_scVI", "X_umap_scANVI", "X_umap"]


# ── Per-cell normalisation ─────────────────────────────────────────────────────
def lognorm(X, target_sum=1e4):
    """normalize_total + log1p on a sparse or dense matrix (cells × genes)."""
    if sp.issparse(X):
        X = X.astype(np.float32)
        totals = np.asarray(X.sum(axis=1)).ravel()
        totals[totals == 0] = 1
        # Divide each row by its total
        diag = sp.diags(1.0 / totals)
        X = diag.dot(X) * target_sum
        X.data = np.log1p(X.data)
    else:
        X = X.astype(np.float32)
        totals = X.sum(axis=1, keepdims=True)
        totals[totals == 0] = 1
        X = np.log1p(X / totals * target_sum)
    return X


# ── Load raw count matrix from one source h5ad ────────────────────────────────
def load_source_counts(h5ad_path, needed_barcodes):
    """
    Return (AnnData with raw counts in X, original obs_names) restricted to
    needed_barcodes.  var_names are gene symbols.
    """
    a = sc.read_h5ad(h5ad_path)

    # CellxGene: remap var_names to gene symbols
    if "feature_name" in a.var.columns:
        a.var_names = a.var["feature_name"].astype(str).values
        a.var_names_make_unique()

    # Priority: raw.X (CellxGene integer counts) → layers["counts"] → X
    if a.raw is not None:
        raw = a.raw.to_adata()
        if "feature_name" in raw.var.columns:
            raw.var_names = raw.var["feature_name"].astype(str).values
            raw.var_names_make_unique()
        # Rebuild minimal AnnData with raw counts
        bc_mask = raw.obs_names.isin(needed_barcodes)
        a_sub   = raw[bc_mask].copy()
        a_sub.X = sp.csr_matrix(a_sub.X).astype(np.float32)
        a_sub.X.data = np.round(a_sub.X.data)
        return a_sub

    # Subset to needed barcodes
    bc_mask = a.obs_names.isin(needed_barcodes)
    a_sub   = a[bc_mask].copy()

    if "counts" in a_sub.layers:
        X = sp.csr_matrix(a_sub.layers["counts"]).astype(np.float32)
    else:
        X = sp.csr_matrix(a_sub.X).astype(np.float32)

    X.data = np.round(X.data)

    # If X looks log-norm (non-integer, max < 20), back-transform to pseudo-counts
    vals = X.data[X.data != 0][:2000]
    if len(vals) > 0 and not np.all(vals == np.floor(vals)) and vals.max() < 20:
        print(f"    X appears log-norm — back-transforming to pseudo-counts")
        X.data = np.round(np.expm1(X.data) * 10)

    result = ad.AnnData(
        X   = X,
        obs = a_sub.obs[[]].copy(),   # keep index only
        var = a_sub.var[[]].copy(),
    )
    result.var_names = a_sub.var_names
    return result


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    # ── Step 1: Load obs + obsm from integrated object ────────────────────────
    print("Loading integrated object obs + obsm ...")
    adata_int = sc.read_h5ad(ANNOT_H5AD, backed="r")
    obs_int   = adata_int.obs.copy()

    # Collect UMAP embeddings
    obsm_int = {}
    for key in OBSM_KEYS:
        if key in adata_int.obsm:
            obsm_int[key] = np.array(adata_int.obsm[key])

    adata_int.file.close()

    # Filter to experimental cells
    exp_mask = obs_int["source"].isin(EXPERIMENTAL_SOURCES)
    obs_exp  = obs_int[exp_mask].copy()
    print(f"  Experimental cells: {len(obs_exp):,} / {len(obs_int):,}")

    # Strip batch suffix to recover original barcodes
    obs_exp["_orig_bc"] = obs_exp.index.str.rsplit("-", n=1).str[0]

    # ── Step 2: Load full-gene counts per source ──────────────────────────────
    print("\nLoading full-gene counts from original h5ad files ...")
    source_adatas = []

    for source in sorted(obs_exp["source"].unique()):
        h5ad_path = os.path.join(DATA_DIR, f"{source}.h5ad")
        if not os.path.exists(h5ad_path):
            print(f"  WARNING: {h5ad_path} not found — skipped")
            continue

        src_obs  = obs_exp[obs_exp["source"] == source]
        needed   = src_obs["_orig_bc"].values
        print(f"  {source}: {len(needed):,} cells ...")

        a_src = load_source_counts(h5ad_path, needed)
        if a_src is None or a_src.n_obs == 0:
            print(f"    No matching barcodes — skipped")
            continue

        # Map original barcodes → integrated barcodes
        orig_to_int = dict(zip(src_obs["_orig_bc"], src_obs.index))
        new_index   = [orig_to_int.get(bc, bc) for bc in a_src.obs_names]
        a_src.obs_names = pd.Index(new_index)

        print(f"    {a_src.n_obs:,} cells × {a_src.n_vars:,} genes")
        source_adatas.append(a_src)
        gc.collect()

    if not source_adatas:
        sys.exit("ERROR: No source files loaded.")

    # ── Step 3: Concatenate with outer join ───────────────────────────────────
    print("\nConcatenating across sources (outer join) ...")
    adata = ad.concat(
        source_adatas,
        join         = "outer",
        merge        = "first",
        index_unique = None,    # barcodes are already unique (integrated format)
    )
    del source_adatas
    gc.collect()

    # Fill NaN from outer join with 0
    adata.X = sp.csr_matrix(adata.X)
    adata.X.data = np.nan_to_num(adata.X.data, nan=0.0)
    print(f"  Concatenated: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # ── Step 4: Reorder to match obs_exp index ────────────────────────────────
    shared = adata.obs_names.intersection(obs_exp.index)
    adata  = adata[shared].copy()
    print(f"  After barcode alignment: {adata.n_obs:,} cells")

    # ── Step 5: Store raw counts + log-norm ──────────────────────────────────
    adata.layers["counts"] = sp.csr_matrix(adata.X).astype(np.float32)
    print("  Computing log-normalisation ...")
    adata.layers["lognorm"] = lognorm(adata.layers["counts"])
    adata.X = adata.layers["lognorm"].copy()

    # ── Step 6: Transfer metadata and embeddings ──────────────────────────────
    print("  Transferring metadata from integrated object ...")
    adata.obs = obs_exp.loc[adata.obs_names].drop(columns=["_orig_bc"])

    for key, emb in obsm_int.items():
        emb_df = pd.DataFrame(
            emb,
            index=obs_int.index,
        )
        shared_emb = emb_df.loc[adata.obs_names]
        adata.obsm[key] = shared_emb.values
        print(f"    Transferred obsm: {key}")

    # ── Step 7: Save ──────────────────────────────────────────────────────────
    print(f"\nSaving → {OUT_PATH}")
    adata.write_h5ad(OUT_PATH)

    print("\n===== Full-gene object summary =====")
    print(f"  Cells  : {adata.n_obs:,}")
    print(f"  Genes  : {adata.n_vars:,}")
    print(f"  Layers : {list(adata.layers.keys())}")
    print(f"  obsm   : {list(adata.obsm.keys())}")
    print(f"\n  Cells per condition:")
    print(adata.obs["condition"].value_counts().to_string())


if __name__ == "__main__":
    main()
