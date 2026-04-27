#!/usr/bin/env python3
"""
plot_gene_umap.py — Plot UMAP colored by a single gene's expression,
even if the gene was not among the 3,000 HVGs used for integration.

Raw counts are fetched from the original per-sample h5ad files, normalised
per cell (normalize_total to 10k + log1p), then joined onto the integrated
UMAP coordinates by cell barcode.  One panel is drawn per condition.

Usage
-----
  python plot_gene_umap.py --gene Upk2
  python plot_gene_umap.py --gene Lrp2 --conditions Healthy UUO rUUO
  python plot_gene_umap.py --gene Nphs1 --sources KidneyHealthy1 KidneyUUO1
  python plot_gene_umap.py --gene Aqp2 --umap_key X_umap_scANVI --cell_type_key scanvi_MKA_label
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR   = "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
DATA_DIR   = os.path.join(BASE_DIR, "RenalUrothelium")
SCRIPT_DIR = os.path.join(BASE_DIR, "RenalUrotheliumScripts")
OUT_DIR    = os.path.join(SCRIPT_DIR, "output")
PLOT_DIR   = os.path.join(OUT_DIR, "plots")
ANNOT_H5AD = os.path.join(OUT_DIR, "RenalUrothelium_integrated_annotated.h5ad")

os.makedirs(PLOT_DIR, exist_ok=True)


# ── Argument parsing ───────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gene", required=True,
                   help="Gene symbol to plot (e.g. Upk2, Lrp2)")
    p.add_argument("--conditions", nargs="*", default=None,
                   help="Conditions to include. Default: all. E.g. --conditions Healthy UUO")
    p.add_argument("--sources", nargs="*", default=None,
                   help="Source sample IDs to include. Default: all. "
                        "E.g. --sources KidneyHealthy1 KidneyUUO1")
    p.add_argument("--umap_key", default="X_umap_scVI",
                   help="obsm key for UMAP coordinates (default: X_umap_scVI)")
    p.add_argument("--cell_type_key", default="scanvi_Lake_label",
                   help="obs column for cell-type labels (default: scanvi_Lake_label)")
    p.add_argument("--out", default=None,
                   help="Output PDF path. Default: output/plots/umap_<gene>.pdf")
    return p.parse_args()


# ── Per-cell normalisation (same as 01_preprocess.py) ─────────────────────────
def lognorm(counts_vec, totals_vec, target_sum=1e4):
    """normalize_total + log1p, per cell."""
    totals_vec = totals_vec.copy().astype(float)
    totals_vec[totals_vec == 0] = 1          # avoid division by zero
    return np.log1p(counts_vec / totals_vec * target_sum)


# ── Extract one gene from a single h5ad file ──────────────────────────────────
def extract_gene(h5ad_path, gene, needed_barcodes):
    """
    Return a pd.Series(index=original_barcode, values=log-norm expression)
    for `gene` from `h5ad_path`, restricted to `needed_barcodes`.
    Returns None if gene not found.
    """
    a = sc.read_h5ad(h5ad_path)

    # ── Remap var_names to gene symbols if CellxGene format ───────────────────
    if gene not in a.var_names and "feature_name" in a.var.columns:
        a.var_names = a.var["feature_name"].astype(str).values
        a.var_names_make_unique()

    # ── Try raw.X first (CellxGene stores raw integer counts there) ───────────
    if a.raw is not None:
        raw = a.raw.to_adata()
        if "feature_name" in raw.var.columns:
            raw.var_names = raw.var["feature_name"].astype(str).values
            raw.var_names_make_unique()
        if gene in raw.var_names:
            idx      = list(raw.var_names).index(gene)
            col      = raw.X[:, idx]
            counts   = np.asarray(col.todense()).ravel() if sp.issparse(col) else np.asarray(col).ravel()
            X_full   = raw.X
            totals   = np.asarray(X_full.sum(axis=1)).ravel() if sp.issparse(X_full) else X_full.sum(axis=1)
            expr     = lognorm(counts, totals)
            series   = pd.Series(expr, index=a.obs_names, name=gene)
            return series.reindex(needed_barcodes).dropna()

    # ── Otherwise use layers["counts"] or X ───────────────────────────────────
    if gene not in a.var_names:
        print(f"    Gene '{gene}' not found in {os.path.basename(h5ad_path)} — skipped")
        return None

    idx = list(a.var_names).index(gene)

    if "counts" in a.layers:
        src  = a.layers["counts"]
        full = a.layers["counts"]
    else:
        src  = a.X
        full = a.X

    col    = src[:, idx]
    counts = np.asarray(col.todense()).ravel() if sp.issparse(col) else np.asarray(col).ravel()
    totals = np.asarray(full.sum(axis=1)).ravel() if sp.issparse(full) else full.sum(axis=1)

    # If X is already log-norm (float, max < 20) use it directly without re-normalising
    vals = counts[counts != 0][:1000]
    if len(vals) > 0 and not np.all(vals == np.floor(vals)) and vals.max() < 20:
        print(f"    {os.path.basename(h5ad_path)}: X appears log-normalised — using directly")
        series = pd.Series(counts, index=a.obs_names, name=gene)
    else:
        expr   = lognorm(counts, totals)
        series = pd.Series(expr, index=a.obs_names, name=gene)

    return series.reindex(needed_barcodes).dropna()


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    args     = parse_args()
    gene     = args.gene
    out_path = args.out or os.path.join(PLOT_DIR, f"umap_{gene}.pdf")

    # ── Load metadata + UMAP from integrated object (backed=r skips loading X) ─
    print(f"Loading integrated object (obs + obsm only) ...")
    adata = sc.read_h5ad(ANNOT_H5AD, backed="r")

    umap_key = args.umap_key
    if umap_key not in adata.obsm:
        fallback = [k for k in adata.obsm if "umap" in k.lower()]
        if not fallback:
            sys.exit(f"ERROR: '{umap_key}' not found in obsm. Available: {list(adata.obsm)}")
        umap_key = fallback[0]
        print(f"  Warning: using fallback UMAP key '{umap_key}'")

    obs = adata.obs.copy()
    umap = pd.DataFrame(
        np.array(adata.obsm[umap_key]),
        index=adata.obs_names,
        columns=["UMAP1", "UMAP2"],
    )
    obs = obs.join(umap)
    adata.file.close()

    print(f"  Total cells: {len(obs):,}")

    # ── Optional filters ───────────────────────────────────────────────────────
    if args.conditions:
        obs = obs[obs["condition"].isin(args.conditions)]
        print(f"  After condition filter {args.conditions}: {len(obs):,} cells")
    if args.sources:
        obs = obs[obs["source"].isin(args.sources)]
        print(f"  After source filter {args.sources}: {len(obs):,} cells")

    # ── Recover original barcodes (strip batch suffix added by ad.concat) ──────
    # Integrated barcode format: <original_barcode>-<batch_index>
    obs["_orig_bc"] = obs.index.str.rsplit("-", n=1).str[0]

    # ── Fetch gene expression per source ──────────────────────────────────────
    print(f"\nFetching '{gene}' from original h5ad files ...")
    expr_parts = []

    for source in sorted(obs["source"].unique()):
        h5ad_path = os.path.join(DATA_DIR, f"{source}.h5ad")
        if not os.path.exists(h5ad_path):
            print(f"  WARNING: {h5ad_path} not found — skipped")
            continue

        src_rows       = obs[obs["source"] == source]
        needed_barcodes = src_rows["_orig_bc"].values
        print(f"  {source}: {len(needed_barcodes):,} cells ...")

        series = extract_gene(h5ad_path, gene, needed_barcodes)
        if series is None or series.empty:
            continue

        # Map original barcode → integrated index
        bc_map  = src_rows[["_orig_bc"]].reset_index()   # cols: index (integrated), _orig_bc
        bc_map  = bc_map.join(series.rename(gene), on="_orig_bc")
        bc_map  = bc_map.set_index("index")
        expr_parts.append(bc_map[gene])

    if not expr_parts:
        sys.exit(f"ERROR: Gene '{gene}' was not found in any source h5ad file.")

    expr_all   = pd.concat(expr_parts)
    obs[gene]  = expr_all
    obs[gene]  = obs[gene].fillna(0.0)

    n_nonzero  = (obs[gene] > 0).sum()
    print(f"\n'{gene}': {n_nonzero:,} / {len(obs):,} cells non-zero "
          f"({100 * n_nonzero / len(obs):.1f}%)")
    print(f"  Expression range: {obs[gene].min():.3f} – {obs[gene].max():.3f}")

    # ── Plot: one panel per condition ─────────────────────────────────────────
    conditions = sorted(obs["condition"].unique())
    ncols      = len(conditions)
    fig, axes  = plt.subplots(1, ncols, figsize=(4.5 * ncols, 4.2), squeeze=False)

    vmax = obs[gene].quantile(0.99) or obs[gene].max()
    vmin = 0.0

    for ax, cond in zip(axes[0], conditions):
        sub = obs[obs["condition"] == cond]

        # Grey background = all cells
        ax.scatter(obs["UMAP1"], obs["UMAP2"],
                   s=0.2, c="lightgrey", alpha=0.4, linewidths=0, rasterized=True)

        # Condition cells coloured by expression
        sc_kw = dict(c=sub[gene], cmap="Reds", vmin=vmin, vmax=vmax,
                     s=0.5, alpha=0.8, linewidths=0, rasterized=True)
        sc_obj = ax.scatter(sub["UMAP1"], sub["UMAP2"], **sc_kw)

        ax.set_title(f"{cond}\n(n={len(sub):,})", fontsize=10)
        ax.axis("off")

    cbar = plt.colorbar(sc_obj, ax=axes[0, -1], shrink=0.7, pad=0.02)
    cbar.set_label(f"{gene}\n(log-norm)", fontsize=9)
    fig.suptitle(f"{gene} expression by condition", fontsize=13, y=1.01)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"\nSaved → {out_path}")


if __name__ == "__main__":
    main()
