#!/usr/bin/env python3
"""
08_pseudobulk.py — Build a pseudobulk count matrix for DEG analysis.

For each experimental sample, raw integer counts are fetched from the
original per-sample h5ad files using cell barcodes from the integrated
object.  Counts are summed per (sample × cell_type) group (pseudobulk).

Only the 14 experimental samples are used — reference atlases (Lake,
MKA, Chen, UUOProjectObject) are excluded because they were included
for integration quality but are not part of the experimental design.

Outputs (in output/pseudobulk/)
--------------------------------
  counts.csv   — genes × pseudobulk-samples (integer sums)
  meta.csv     — pseudobulk-sample metadata: sample_id, condition,
                 cell_type, n_cells

Usage
-----
  python 08_pseudobulk.py
  python 08_pseudobulk.py --cell_type_key scanvi_MKA_label
  python 08_pseudobulk.py --cell_types PT DCT TAL   # subset cell types
  python 08_pseudobulk.py --min_cells 20            # min cells per pseudobulk
"""

import argparse
import os
import gc
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR   = "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
DATA_DIR   = os.path.join(BASE_DIR, "RenalUrothelium")
SCRIPT_DIR = os.path.join(BASE_DIR, "RenalUrotheliumScripts")
OUT_DIR    = os.path.join(SCRIPT_DIR, "output")
ANNOT_H5AD = os.path.join(OUT_DIR, "RenalUrothelium_integrated_annotated.h5ad")
PB_DIR     = os.path.join(OUT_DIR, "pseudobulk")

os.makedirs(PB_DIR, exist_ok=True)

# ── Experimental samples only (exclude reference atlases) ─────────────────────
EXPERIMENTAL_SOURCES = {
    "KidneyHealthy1", "KidneyHealthy2", "KidneyHealthy3",
    "KidneyHealthy4", "KidneyHealthy5",
    "KidneyTET2UUO",
    "KidneyUUO1",  "KidneyUUO2",  "KidneyUUO3",  "KidneyUUO4",
    "KidneyUUO5",  "KidneyUUO6",  "KidneyUUO7",  "KidneyUUO8",
    "KidneyrUUO1",
}


# ── Argument parsing ───────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cell_type_key", default="scanvi_Lake_label",
                   help="obs column for cell-type labels (default: scanvi_Lake_label). "
                        "Options: scanvi_Lake_label, scanvi_MKA_label, leiden_scVI")
    p.add_argument("--cell_types", nargs="*", default=None,
                   help="Cell types to include. Default: all. "
                        "E.g. --cell_types PT DCT TAL PC IC")
    p.add_argument("--min_cells", type=int, default=10,
                   help="Minimum cells per pseudobulk sample (default: 10). "
                        "Groups below this threshold are dropped.")
    p.add_argument("--out_dir", default=PB_DIR,
                   help=f"Output directory (default: {PB_DIR})")
    return p.parse_args()


# ── Helper: extract raw counts for one gene set from a source file ─────────────
def get_raw_count_matrix(h5ad_path, needed_barcodes):
    """
    Load the source h5ad and return (sparse count matrix, gene names, obs_names)
    restricted to needed_barcodes.  Returns raw integer counts.
    """
    a = sc.read_h5ad(h5ad_path)

    # CellxGene format: remap var_names to gene symbols
    if "feature_name" in a.var.columns:
        a.var_names = a.var["feature_name"].astype(str).values
        a.var_names_make_unique()

    # Priority: raw.X (CellxGene integer counts) → layers["counts"] → X
    if a.raw is not None:
        raw = a.raw.to_adata()
        if "feature_name" in raw.var.columns:
            raw.var_names = raw.var["feature_name"].astype(str).values
            raw.var_names_make_unique()
        X      = raw.X
        genes  = list(raw.var_names)
        obs_bc = list(a.obs_names)       # raw shares obs with a
    elif "counts" in a.layers:
        X      = a.layers["counts"]
        genes  = list(a.var_names)
        obs_bc = list(a.obs_names)
    else:
        X      = a.X
        genes  = list(a.var_names)
        obs_bc = list(a.obs_names)

    # Subset to needed barcodes
    bc_index   = {bc: i for i, bc in enumerate(obs_bc)}
    row_idx    = [bc_index[bc] for bc in needed_barcodes if bc in bc_index]
    found_bcs  = [needed_barcodes[i] for i, bc in enumerate(needed_barcodes)
                  if bc in bc_index]

    if not row_idx:
        return None, None, None

    X_sub = X[row_idx, :]
    if not sp.issparse(X_sub):
        X_sub = sp.csr_matrix(X_sub)
    else:
        X_sub = sp.csr_matrix(X_sub)

    # Round to integers (some sources store float counts)
    X_sub.data = np.round(X_sub.data).astype(np.int32)

    return X_sub, genes, found_bcs


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # ── Load obs from integrated object ───────────────────────────────────────
    print("Loading integrated object obs ...")
    adata = sc.read_h5ad(ANNOT_H5AD, backed="r")
    obs   = adata.obs[[
        "source", "sample_id", "condition",
        args.cell_type_key,
    ]].copy()
    obs.index.name = "barcode"
    adata.file.close()

    obs = obs.rename(columns={args.cell_type_key: "cell_type"})
    obs["cell_type"] = obs["cell_type"].astype(str)

    # ── Filter to experimental samples only ───────────────────────────────────
    obs = obs[obs["source"].isin(EXPERIMENTAL_SOURCES)]
    print(f"  Experimental cells: {len(obs):,}")

    # ── Optional cell-type filter ─────────────────────────────────────────────
    if args.cell_types:
        obs = obs[obs["cell_type"].isin(args.cell_types)]
        print(f"  After cell-type filter {args.cell_types}: {len(obs):,} cells")

    # Strip batch suffix → original barcode
    obs["_orig_bc"] = obs.index.str.rsplit("-", n=1).str[0]

    # ── Load raw counts per source and accumulate ──────────────────────────────
    # We build a list of (cell_metadata, count_matrix, gene_list) per source,
    # then align all to a common gene set and pseudobulk-aggregate.
    print("\nLoading raw counts from original h5ad files ...")

    per_source = []   # list of (sub_obs, X_csr, genes)

    for source in sorted(obs["source"].unique()):
        h5ad_path = os.path.join(DATA_DIR, f"{source}.h5ad")
        if not os.path.exists(h5ad_path):
            print(f"  WARNING: {h5ad_path} not found — skipped")
            continue

        src_obs = obs[obs["source"] == source]
        needed  = src_obs["_orig_bc"].values
        print(f"  {source}: {len(needed):,} cells ...")

        X_sub, genes, found_bcs = get_raw_count_matrix(h5ad_path, needed)
        if X_sub is None:
            print(f"    No matching barcodes — skipped")
            continue

        # Rebuild sub_obs aligned to found_bcs
        orig_to_integrated = src_obs.reset_index().set_index("_orig_bc")["barcode"]
        integrated_bcs     = [orig_to_integrated[bc] for bc in found_bcs
                               if bc in orig_to_integrated.index]
        sub_obs = src_obs.loc[integrated_bcs].copy()
        sub_obs["_found_bc"] = found_bcs[:len(integrated_bcs)]

        per_source.append((sub_obs, X_sub, genes))
        gc.collect()

    if not per_source:
        raise RuntimeError("No source files could be loaded.")

    # ── Build union gene set ───────────────────────────────────────────────────
    print("\nBuilding union gene set ...")
    all_genes = []
    for _, _, genes in per_source:
        all_genes.extend(genes)
    union_genes = list(dict.fromkeys(all_genes))   # preserve order, deduplicate
    gene_index  = {g: i for i, g in enumerate(union_genes)}
    n_genes     = len(union_genes)
    print(f"  Union gene set: {n_genes:,} genes")

    # ── Pseudobulk aggregation ─────────────────────────────────────────────────
    print("\nAggregating pseudobulk counts ...")

    # Group key: sample_id × cell_type
    all_obs = pd.concat([so for so, _, _ in per_source])
    groups  = all_obs.groupby(["sample_id", "cell_type"]).size()
    groups  = groups[groups >= args.min_cells]
    print(f"  Pseudobulk groups (≥{args.min_cells} cells): {len(groups)}")

    pb_counts = {}   # key → dense count vector (n_genes,)
    pb_meta   = []

    for source_idx, (sub_obs, X_sub, genes) in enumerate(per_source):
        # Map this source's gene indices into the union gene index
        col_map = np.array([gene_index[g] for g in genes], dtype=np.int32)

        for (sample_id, cell_type), grp in sub_obs.groupby(["sample_id", "cell_type"]):
            key = f"{sample_id}__{cell_type}"
            if (sample_id, cell_type) not in groups.index:
                continue   # below min_cells threshold

            # Row indices in X_sub for this group
            group_rows = [i for i, bc in enumerate(sub_obs.index)
                          if bc in grp.index]
            if not group_rows:
                continue

            X_grp = X_sub[group_rows, :]
            # Sum counts across cells
            summed = np.asarray(X_grp.sum(axis=0)).ravel()   # (n_source_genes,)

            if key not in pb_counts:
                pb_counts[key] = np.zeros(n_genes, dtype=np.int64)
            pb_counts[key][col_map] += summed.astype(np.int64)

            if key not in {m["pseudobulk_id"] for m in pb_meta}:
                cond = grp["condition"].iloc[0]
                pb_meta.append({
                    "pseudobulk_id": key,
                    "sample_id":     sample_id,
                    "condition":     cond,
                    "cell_type":     cell_type,
                    "n_cells":       len(grp),
                })

    # ── Save ──────────────────────────────────────────────────────────────────
    counts_path = os.path.join(args.out_dir, "counts.csv")
    meta_path   = os.path.join(args.out_dir, "meta.csv")

    print(f"\nSaving {n_genes:,} genes × {len(pb_counts)} pseudobulk samples ...")

    counts_df = pd.DataFrame(pb_counts, index=union_genes)
    counts_df.index.name = "gene"
    counts_df.to_csv(counts_path)
    print(f"  Counts → {counts_path}")

    meta_df = pd.DataFrame(pb_meta).set_index("pseudobulk_id")
    meta_df.to_csv(meta_path)
    print(f"  Meta   → {meta_path}")

    print("\n===== Pseudobulk summary =====")
    print(meta_df["condition"].value_counts().to_string())
    print(f"\nCell types:")
    print(meta_df["cell_type"].value_counts().head(20).to_string())


if __name__ == "__main__":
    main()
