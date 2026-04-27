#!/usr/bin/env python3
"""
plot_gene_umap.py — Plot UMAP colored by a single gene's expression,
corrected for sample and technology batch effects where possible.

Strategy
--------
  HVG (in the 3,000 used for integration)
      → scVI get_normalized_expression(): decoder-corrected, removes
        sample and technology batch effects.  Requires loading the
        scVI model (~few GB RAM).

  Non-HVG
      → Raw counts fetched from original per-sample h5ad files,
        normalised per cell (normalize_total 10k + log1p).
        Not batch-corrected at the expression level, but the UMAP
        positions are already scVI-corrected.

One panel is drawn per condition.

Usage
-----
  python plot_gene_umap.py --gene Upk2
  python plot_gene_umap.py --gene Lrp2 --conditions Healthy UUO rUUO
  python plot_gene_umap.py --gene Nphs1 --sources KidneyHealthy1 KidneyUUO1
  python plot_gene_umap.py --gene Aqp2 --umap_key X_umap_scANVI
  python plot_gene_umap.py --gene Lrp2 --no_scvi   # force raw log-norm even for HVGs
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
BASE_DIR    = "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
DATA_DIR    = os.path.join(BASE_DIR, "RenalUrothelium")
SCRIPT_DIR  = os.path.join(BASE_DIR, "RenalUrotheliumScripts")
OUT_DIR     = os.path.join(SCRIPT_DIR, "output")
PLOT_DIR    = os.path.join(OUT_DIR, "plots")
ANNOT_H5AD  = os.path.join(OUT_DIR, "RenalUrothelium_integrated_annotated.h5ad")
SCVI_MODEL  = os.path.join(OUT_DIR, "scvi_model")

os.makedirs(PLOT_DIR, exist_ok=True)


# ── Argument parsing ───────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gene", required=True,
                   help="Gene symbol to plot (e.g. Upk2, Lrp2)")
    p.add_argument("--conditions", nargs="*", default=None,
                   help="Conditions to include. Default: all. "
                        "E.g. --conditions Healthy UUO rUUO")
    p.add_argument("--sources", nargs="*", default=None,
                   help="Source sample IDs to include. Default: all. "
                        "E.g. --sources KidneyHealthy1 KidneyUUO1")
    p.add_argument("--umap_key", default="X_umap_scVI",
                   help="obsm key for UMAP coordinates (default: X_umap_scVI)")
    p.add_argument("--no_scvi", action="store_true",
                   help="Force raw log-norm even for HVGs (skip scVI model loading)")
    p.add_argument("--out", default=None,
                   help="Output PDF path. Default: output/plots/umap_<gene>.pdf")
    return p.parse_args()


# ── scVI normalised expression (HVGs only, batch-corrected) ───────────────────
def get_scvi_expression(gene, obs_index):
    """
    Load scVI model and return decoder-normalised expression for `gene`.
    Accounts for sample and technology batch effects.
    Returns pd.Series(index=cell_barcode) or None on failure.
    """
    try:
        import scvi
    except ImportError:
        print("  scvi-tools not available — falling back to raw log-norm")
        return None

    if not os.path.isdir(SCVI_MODEL):
        print(f"  scVI model not found at {SCVI_MODEL} — falling back to raw log-norm")
        return None

    print(f"  Loading scVI model from {SCVI_MODEL} ...")
    # Load full adata (counts layer needed for encoding)
    adata = sc.read_h5ad(ANNOT_H5AD)

    if gene not in adata.var_names:
        print(f"  '{gene}' not in HVGs — this branch should not be reached")
        return None

    try:
        model = scvi.model.SCVI.load(SCVI_MODEL, adata=adata)
        print(f"  Computing scVI normalised expression for '{gene}' ...")
        expr_df = model.get_normalized_expression(
            gene_list=[gene],
            library_size="latent",   # uses model's latent library size estimate
            n_samples=1,
        )
        series = expr_df[gene]
        series.index = adata.obs_names
        # log1p so scale matches the raw log-norm fallback
        series = np.log1p(series * 1e4)
        return series.reindex(obs_index)
    except Exception as e:
        print(f"  scVI normalisation failed ({e}) — falling back to raw log-norm")
        return None


# ── Per-cell normalisation (same as 01_preprocess.py) ─────────────────────────
def lognorm(counts_vec, totals_vec, target_sum=1e4):
    totals_vec = totals_vec.copy().astype(float)
    totals_vec[totals_vec == 0] = 1
    return np.log1p(counts_vec / totals_vec * target_sum)


# ── Extract one gene from a single original h5ad file ─────────────────────────
def extract_gene_from_source(h5ad_path, gene, needed_barcodes):
    """
    Return pd.Series(index=original_barcode, values=log-norm expression).
    Returns None if gene not found in the file.
    """
    a = sc.read_h5ad(h5ad_path)

    # Remap var_names to gene symbols for CellxGene-format files
    if gene not in a.var_names and "feature_name" in a.var.columns:
        a.var_names = a.var["feature_name"].astype(str).values
        a.var_names_make_unique()

    # CellxGene: raw.X has integer counts over the full gene set
    if a.raw is not None:
        raw = a.raw.to_adata()
        if "feature_name" in raw.var.columns:
            raw.var_names = raw.var["feature_name"].astype(str).values
            raw.var_names_make_unique()
        if gene in raw.var_names:
            idx    = list(raw.var_names).index(gene)
            col    = raw.X[:, idx]
            counts = np.asarray(col.todense()).ravel() if sp.issparse(col) else np.asarray(col).ravel()
            totals = np.asarray(raw.X.sum(axis=1)).ravel() if sp.issparse(raw.X) else raw.X.sum(axis=1)
            series = pd.Series(lognorm(counts, totals), index=a.obs_names, name=gene)
            return series.reindex(needed_barcodes).dropna()

    if gene not in a.var_names:
        print(f"    '{gene}' not found in {os.path.basename(h5ad_path)} — skipped")
        return None

    idx = list(a.var_names).index(gene)
    src = a.layers["counts"] if "counts" in a.layers else a.X
    col = src[:, idx]
    counts = np.asarray(col.todense()).ravel() if sp.issparse(col) else np.asarray(col).ravel()
    totals = np.asarray(src.sum(axis=1)).ravel() if sp.issparse(src) else src.sum(axis=1)

    # If X is already log-norm (non-integer floats with max < 20), use directly
    vals = counts[counts != 0][:1000]
    if len(vals) > 0 and not np.all(vals == np.floor(vals)) and vals.max() < 20:
        print(f"    {os.path.basename(h5ad_path)}: X appears log-normalised — using directly")
        series = pd.Series(counts, index=a.obs_names, name=gene)
    else:
        series = pd.Series(lognorm(counts, totals), index=a.obs_names, name=gene)

    return series.reindex(needed_barcodes).dropna()


# ── Fetch expression via raw log-norm from original files ──────────────────────
def get_raw_lognorm_expression(gene, obs):
    """
    For each source in obs, open the original h5ad and extract log-norm
    expression for `gene`.  Returns obs with a new column for the gene.
    """
    obs = obs.copy()
    obs["_orig_bc"] = obs.index.str.rsplit("-", n=1).str[0]
    expr_parts = []

    for source in sorted(obs["source"].unique()):
        h5ad_path = os.path.join(DATA_DIR, f"{source}.h5ad")
        if not os.path.exists(h5ad_path):
            print(f"  WARNING: {h5ad_path} not found — skipped")
            continue

        src_rows = obs[obs["source"] == source]
        print(f"  {source}: {len(src_rows):,} cells ...")

        series = extract_gene_from_source(h5ad_path, gene, src_rows["_orig_bc"].values)
        if series is None or series.empty:
            continue

        bc_map = src_rows[["_orig_bc"]].reset_index()
        bc_map = bc_map.join(series.rename(gene), on="_orig_bc")
        bc_map = bc_map.set_index("index")
        expr_parts.append(bc_map[gene])

    if not expr_parts:
        sys.exit(f"ERROR: Gene '{gene}' was not found in any source h5ad file.")

    expr_all  = pd.concat(expr_parts)
    obs[gene] = expr_all
    obs[gene] = obs[gene].fillna(0.0)
    return obs


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    args     = parse_args()
    gene     = args.gene
    out_path = args.out or os.path.join(PLOT_DIR, f"umap_{gene}.pdf")

    # ── Load metadata + UMAP (backed mode avoids loading 31 GB X) ─────────────
    print("Loading integrated object metadata ...")
    adata = sc.read_h5ad(ANNOT_H5AD, backed="r")

    umap_key = args.umap_key
    if umap_key not in adata.obsm:
        fallback = [k for k in adata.obsm if "umap" in k.lower()]
        if not fallback:
            sys.exit(f"ERROR: '{umap_key}' not found. Available: {list(adata.obsm)}")
        umap_key = fallback[0]
        print(f"  Warning: using fallback UMAP key '{umap_key}'")

    hvg_genes = list(adata.var_names)           # the 3,000 HVGs
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

    # ── Choose expression strategy ─────────────────────────────────────────────
    is_hvg     = gene in hvg_genes
    use_scvi   = is_hvg and not args.no_scvi and os.path.isdir(SCVI_MODEL)
    expr_label = gene  # label for colorbar

    if use_scvi:
        print(f"\n'{gene}' is an HVG — using scVI normalised expression (batch-corrected) ...")
        scvi_expr = get_scvi_expression(gene, obs.index)
        if scvi_expr is not None:
            obs[gene]  = scvi_expr.reindex(obs.index).fillna(0.0)
            expr_label = f"{gene}\n(scVI-normalised)"
        else:
            use_scvi = False   # fallback triggered inside get_scvi_expression

    if not use_scvi:
        if is_hvg:
            print(f"\n'{gene}' is an HVG but using raw log-norm (scVI unavailable or --no_scvi set)")
        else:
            print(f"\n'{gene}' is NOT an HVG — fetching from original h5ad files ...")
        obs        = get_raw_lognorm_expression(gene, obs)
        expr_label = f"{gene}\n(log-norm, not batch-corrected)"

    n_nonzero = (obs[gene] > 0).sum()
    print(f"\n'{gene}': {n_nonzero:,} / {len(obs):,} cells non-zero "
          f"({100 * n_nonzero / len(obs):.1f}%)")
    print(f"  Range: {obs[gene].min():.3f} – {obs[gene].max():.3f}")

    # ── Plot: one panel per condition ──────────────────────────────────────────
    conditions = sorted(obs["condition"].unique())
    ncols      = len(conditions)
    fig, axes  = plt.subplots(1, ncols, figsize=(4.5 * ncols, 4.2), squeeze=False)

    vmax = obs[gene].quantile(0.99) or obs[gene].max()
    vmin = 0.0

    for ax, cond in zip(axes[0], conditions):
        sub = obs[obs["condition"] == cond]

        # Grey background — all cells for spatial context
        ax.scatter(obs["UMAP1"], obs["UMAP2"],
                   s=0.2, c="lightgrey", alpha=0.4, linewidths=0, rasterized=True)
        sc_obj = ax.scatter(sub["UMAP1"], sub["UMAP2"],
                            c=sub[gene], cmap="Reds",
                            vmin=vmin, vmax=vmax,
                            s=0.5, alpha=0.8, linewidths=0, rasterized=True)
        ax.set_title(f"{cond}\n(n={len(sub):,})", fontsize=10)
        ax.axis("off")

    cbar = plt.colorbar(sc_obj, ax=axes[0, -1], shrink=0.7, pad=0.02)
    cbar.set_label(expr_label, fontsize=9)
    strategy = "scVI batch-corrected" if use_scvi else "raw log-norm"
    fig.suptitle(f"{gene} expression by condition  [{strategy}]", fontsize=12, y=1.01)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"\nSaved → {out_path}")


if __name__ == "__main__":
    main()
