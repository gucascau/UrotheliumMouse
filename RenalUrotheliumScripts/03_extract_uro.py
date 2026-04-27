#!/usr/bin/env python3
"""
03_extract_uro.py — Strictly extract urothelial cells from the
RenalUrothelium scVI-integrated AnnData object.

Filtering logic mirrors the stricter Seurat-style rule:

  Epcam > 0
  AND at least one of Upk1a / Upk1b / Upk2 / Upk3a > 0
  AND KidneyEpiScore1 < 0.1

KidneyEpiScore1 is computed with scanpy.tl.score_genes using kidney epithelial
exclusion markers, similar in spirit to Seurat AddModuleScore.

Input  : output/RenalUrothelium_allcells_scvi.h5ad
Output : output/RenalUrothelium_uro_cells_fullgene_scvi.h5ad
"""

import os

# Keep Scanpy/Numba/Matplotlib caches on writable storage during Slurm runs.
os.environ.setdefault("MPLCONFIGDIR", f"/tmp/matplotlib-{os.getuid()}")
os.environ.setdefault("NUMBA_CACHE_DIR", f"/tmp/numba-{os.getuid()}")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)
os.makedirs(os.environ["NUMBA_CACHE_DIR"], exist_ok=True)

import inspect
import numpy as np
import scanpy as sc
import scipy.sparse as sp

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR = "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
SCRIPT_DIR = os.path.join(BASE_DIR, "RenalUrotheliumScripts")
OUT_DIR = os.path.join(SCRIPT_DIR, "output")

IN_PATH = os.path.join(OUT_DIR, "RenalUrothelium_allcells_scvi.h5ad")
OUT_PATH = os.path.join(OUT_DIR, "RenalUrothelium_uro_cells_fullgene_scvi.h5ad")

# ── Markers ───────────────────────────────────────────────────────────────────
URO_MARKERS_MOUSE = [
    "Krt8", "Krt18", "Krt19", "Epcam",
    "Upk1a", "Upk1b", "Upk2", "Upk3a",
    "Krt20", "Krt5", "Krt14", "Trp63",
    "Foxa1", "Gata3", "Pparg",
]

UPK_MARKERS_MOUSE = ["Upk1a", "Upk1b", "Upk2", "Upk3a"]

KIDNEY_EXCLUDE_MOUSE = [
    "Slc34a1", "Lrp2", "Cubn",      # proximal tubule
    "Umod", "Slc12a1",              # TAL
    "Slc12a3",                      # DCT
    "Aqp2", "Aqp3", "Scnn1g",       # collecting duct principal
    "Atp6v1b1", "Slc4a1", "Foxi1",  # intercalated cells
]

KIDNEY_SCORE_NAME = "KidneyEpiScore1"
KIDNEY_SCORE_MAX = 0.1
RANDOM_STATE = 0


def split_present_missing(adata, genes):
    """Return genes present and absent from adata.var_names, preserving order."""
    present = [g for g in genes if g in adata.var_names]
    missing = [g for g in genes if g not in adata.var_names]
    return present, missing


def expression_for_genes(adata, genes):
    """Extract adata.X columns for genes."""
    idx = [adata.var_names.get_loc(g) for g in genes]
    return adata.X[:, idx]


def positive_any_gene_mask(expr):
    """Cells with expression > 0 in at least one supplied gene."""
    if sp.issparse(expr):
        return np.asarray((expr > 0).sum(axis=1)).ravel() > 0
    return np.asarray(expr > 0).any(axis=1)


def positive_counts_by_gene(expr):
    """Number of cells with expression > 0 for each supplied gene column."""
    if sp.issparse(expr):
        return np.asarray((expr > 0).sum(axis=0)).ravel().astype(int)
    return np.asarray((expr > 0).sum(axis=0)).ravel().astype(int)


def score_genes_like_add_module_score(adata, genes, score_name):
    """Compute a Scanpy module score with stable kwargs across versions."""
    ctrl_size = min(50, max(1, adata.n_vars - len(genes)))
    kwargs = {
        "gene_list": genes,
        "score_name": score_name,
        "ctrl_size": ctrl_size,
        "random_state": RANDOM_STATE,
        "use_raw": False,
    }

    signature = inspect.signature(sc.tl.score_genes)
    supported_kwargs = {k: v for k, v in kwargs.items() if k in signature.parameters}
    sc.tl.score_genes(adata, **supported_kwargs)


def print_gate_count(label, mask, total):
    n = int(mask.sum())
    print(f"  {label:<45s}: {n:>10,} / {total:,} ({100 * n / total:6.2f}%)")


###############################################################################
# Load
###############################################################################

print(f"Loading {IN_PATH} ...")
adata = sc.read_h5ad(IN_PATH)
print(f"  Total: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

###############################################################################
# Marker check
###############################################################################

uro_found, uro_missing = split_present_missing(adata, URO_MARKERS_MOUSE)
upk_found, upk_missing = split_present_missing(adata, UPK_MARKERS_MOUSE)
kidney_found, kidney_missing = split_present_missing(adata, KIDNEY_EXCLUDE_MOUSE)

print(f"\n  Urothelial markers found  ({len(uro_found)}): {uro_found}")
if uro_missing:
    print(f"  Urothelial markers absent ({len(uro_missing)}): {uro_missing}")

print(f"  Uroplakin markers found   ({len(upk_found)}): {upk_found}")
if upk_missing:
    print(f"  Uroplakin markers absent  ({len(upk_missing)}): {upk_missing}")

print(f"  Kidney exclude found      ({len(kidney_found)}): {kidney_found}")
if kidney_missing:
    print(f"  Kidney exclude absent     ({len(kidney_missing)}): {kidney_missing}")

if "Epcam" not in adata.var_names:
    raise ValueError(
        "Strict extraction requires Epcam, but Epcam is absent from this AnnData. "
        "Regenerate the integrated object with Epcam retained, or run this filter "
        "on a full-gene object before HVG restriction."
    )
if not upk_found:
    raise ValueError("Strict extraction requires at least one Upk marker, but none are present.")
if not kidney_found:
    raise ValueError("Strict extraction requires kidney-exclusion markers, but none are present.")

###############################################################################
# Score kidney epithelial exclusion markers
###############################################################################

print("\nScoring kidney epithelial exclusion markers ...")
score_genes_like_add_module_score(adata, kidney_found, KIDNEY_SCORE_NAME)
score = adata.obs[KIDNEY_SCORE_NAME].to_numpy()
print(
    f"  {KIDNEY_SCORE_NAME}: min={np.nanmin(score):.3f}, "
    f"median={np.nanmedian(score):.3f}, max={np.nanmax(score):.3f}"
)

###############################################################################
# Build stringent mask
###############################################################################

epcam_expr = expression_for_genes(adata, ["Epcam"])
upk_expr = expression_for_genes(adata, upk_found)
uro_expr = expression_for_genes(adata, uro_found)

epcam_positive_mask = positive_any_gene_mask(epcam_expr)
upk_positive_mask = positive_any_gene_mask(upk_expr)
kidney_low_mask = score < KIDNEY_SCORE_MAX
strict_uro_mask = epcam_positive_mask & upk_positive_mask & kidney_low_mask

adata.obs["strict_uro_Epcam_positive"] = epcam_positive_mask
adata.obs["strict_uro_Upk_positive"] = upk_positive_mask
adata.obs["strict_uro_low_kidney_epi_score"] = kidney_low_mask
adata.obs["strict_urothelium"] = strict_uro_mask

print("\nStrict urothelial extraction gates:")
print_gate_count("Epcam > 0", epcam_positive_mask, adata.n_obs)
print_gate_count("any Upk1a/Upk1b/Upk2/Upk3a > 0", upk_positive_mask, adata.n_obs)
print_gate_count(f"{KIDNEY_SCORE_NAME} < {KIDNEY_SCORE_MAX}", kidney_low_mask, adata.n_obs)
print_gate_count("ALL strict gates", strict_uro_mask, adata.n_obs)

if int(strict_uro_mask.sum()) == 0:
    raise ValueError("No cells passed the strict urothelial extraction gates.")

###############################################################################
# Marker and metadata summaries
###############################################################################

print("\n  Per-urothelial-marker positive counts:")
for gene, n in zip(uro_found, positive_counts_by_gene(uro_expr)):
    print(f"    {gene:<10s}: {n:>10,}")

print("\n  Per-kidney-exclusion-marker positive counts:")
kidney_expr = expression_for_genes(adata, kidney_found)
for gene, n in zip(kidney_found, positive_counts_by_gene(kidney_expr)):
    print(f"    {gene:<10s}: {n:>10,}")

if "sample_id" in adata.obs.columns:
    print("\n  Strict urothelial cells per sample:")
    print(adata.obs.loc[strict_uro_mask, "sample_id"].value_counts().to_string())

if "condition" in adata.obs.columns:
    print("\n  Strict urothelial cells per condition:")
    print(adata.obs.loc[strict_uro_mask, "condition"].value_counts().to_string())

if "DataSet" in adata.obs.columns:
    print("\n  Strict urothelial cells per DataSet:")
    print(adata.obs.loc[strict_uro_mask, "DataSet"].value_counts().to_string())

if "leiden_scVI" in adata.obs.columns:
    print("\n  Strict urothelial cells per scVI cluster:")
    print(adata.obs.loc[strict_uro_mask, "leiden_scVI"].value_counts().sort_index().to_string())

###############################################################################
# Subset and save
###############################################################################

adata_uro = adata[strict_uro_mask].copy()
print(f"\nSaving {adata_uro.n_obs:,} strict urothelial cells -> {OUT_PATH}")
adata_uro.write_h5ad(OUT_PATH)
print("Done.")
