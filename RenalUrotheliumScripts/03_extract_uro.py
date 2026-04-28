#!/usr/bin/env python3
"""
03_extract_uro.py — Extract urothelial cells from the
RenalUrothelium scVI-integrated AnnData object.

Filtering logic (three-arm OR gate):

  (any of Upk1a/Upk1b/Upk2/Upk3a/Upk3b > 0                    [umbrella arm]
   OR ≥2 of Krt5/Krt14/Trp63 > 0                               [basal arm]
   OR Epcam > 0 AND Krt7 > 0 AND Krt8 > 0 AND any of Foxa1/Gata3 > 0) [intermediate arm]
  AND KidneyEpiScore1 < 0.2

Arm rationale:
  umbrella    — UPK marks differentiated superficial cells
  basal       — ≥2/3 (not all-3) tolerates scRNA-seq dropout in Krt5/Krt14/Trp63
  intermediate — Foxa1/Gata3 are urothelial master TFs absent from all kidney
                 tubular epithelial lineages; combined with Epcam this is highly
                 specific and recovers cells between basal and umbrella layers

KidneyEpiScore1 is computed with scanpy.tl.score_genes using kidney epithelial
exclusion markers, similar in spirit to Seurat AddModuleScore.

Input  : output/RenalUrothelium_allcells_scvi_annotations.h5ad
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
    "Krt8", "Krt18", "Krt19", # epithelial keractin
    "Upk1a", "Upk1b", "Upk2", "Upk3a","Upk3b", # uroplatinks 
    "Krt20", "Krt5", "Krt14", "Trp63", # urothelial subtypes
    "Foxa1", "Gata3", "Pparg",# urothelial transcription factors
]

UPK_MARKERS_MOUSE   = ["Upk1a", "Upk1b", "Upk2", "Upk3a", "Upk3b"]
BASAL_MARKERS_MOUSE = ["Krt5", "Krt14", "Trp63"]   # ≥2/3 required to tolerate dropout
LUMINAL_TF_MOUSE    = ["Foxa1", "Gata3"]            # urothelial master TFs; absent from kidney tubules
BASAL_MIN_POSITIVE  = 2                             # minimum of the 3 basal markers that must be > 0

KIDNEY_EXCLUDE_MOUSE = [
    "Slc34a1", "Lrp2", "Cubn",      # proximal tubule
    "Umod", "Slc12a1",              # TAL
    "Slc12a3",                      # DCT
    "Aqp2", "Aqp3", "Scnn1g",       # collecting duct principal
    "Atp6v1b1", "Slc4a1", "Foxi1",  # intercalated cells
]

KIDNEY_SCORE_NAME = "KidneyEpiScore1"
KIDNEY_SCORE_MAX = 0.2
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


def positive_n_genes_mask(expr, n):
    """Cells with expression > 0 in at least n of the supplied genes."""
    if sp.issparse(expr):
        return np.asarray((expr > 0).sum(axis=1)).ravel() >= n
    return np.asarray((expr > 0).sum(axis=1)).ravel() >= n


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

uro_found, uro_missing       = split_present_missing(adata, URO_MARKERS_MOUSE)
upk_found, upk_missing       = split_present_missing(adata, UPK_MARKERS_MOUSE)
basal_found, basal_missing   = split_present_missing(adata, BASAL_MARKERS_MOUSE)
luminal_found, luminal_missing = split_present_missing(adata, LUMINAL_TF_MOUSE)
kidney_found, kidney_missing = split_present_missing(adata, KIDNEY_EXCLUDE_MOUSE)

print(f"\n  Urothelial markers found  ({len(uro_found)}): {uro_found}")
if uro_missing:
    print(f"  Urothelial markers absent ({len(uro_missing)}): {uro_missing}")
print(f"  Uroplakin markers found   ({len(upk_found)}): {upk_found}")
if upk_missing:
    print(f"  Uroplakin markers absent  ({len(upk_missing)}): {upk_missing}")
print(f"  Basal markers found       ({len(basal_found)}): {basal_found}")
if basal_missing:
    print(f"  Basal markers absent      ({len(basal_missing)}): {basal_missing}")
print(f"  Luminal TF markers found  ({len(luminal_found)}): {luminal_found}")
if luminal_missing:
    print(f"  Luminal TF markers absent ({len(luminal_missing)}): {luminal_missing}")
print(f"  Kidney exclude found      ({len(kidney_found)}): {kidney_found}")
if kidney_missing:
    print(f"  Kidney exclude absent     ({len(kidney_missing)}): {kidney_missing}")

if not upk_found:
    raise ValueError("Extraction requires at least one Upk marker, but none are present.")
if len(basal_found) < BASAL_MIN_POSITIVE:
    raise ValueError(
        f"Basal arm requires at least {BASAL_MIN_POSITIVE} of {BASAL_MARKERS_MOUSE} "
        f"present in the dataset, but only found: {basal_found}"
    )

if "Epcam" not in adata.var_names:
    raise ValueError("Intermediate arm requires Epcam, but it is absent from this AnnData.")
for _g in ("Krt7", "Krt8"):
    if _g not in adata.var_names:
        raise ValueError(f"Intermediate arm requires {_g}, but it is absent from this AnnData.")
if not luminal_found:
    raise ValueError("Intermediate arm requires at least one of Foxa1/Gata3, but none are present.")
if not kidney_found:
    raise ValueError("Extraction requires kidney-exclusion markers, but none are present.")

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

upk_expr     = expression_for_genes(adata, upk_found)
basal_expr   = expression_for_genes(adata, basal_found)
epcam_expr   = expression_for_genes(adata, ["Epcam"])
luminal_expr = expression_for_genes(adata, luminal_found)
uro_expr     = expression_for_genes(adata, uro_found)

# arm 1: umbrella — any UPK > 0
upk_mask          = positive_any_gene_mask(upk_expr)
# arm 2: basal — ≥2 of Krt5/Krt14/Trp63 > 0
basal_mask        = positive_n_genes_mask(basal_expr, BASAL_MIN_POSITIVE)
# arm 3: intermediate — Epcam > 0 AND Krt7 > 0 AND Krt8 > 0 AND any of Foxa1/Gata3 > 0
epcam_mask        = positive_any_gene_mask(epcam_expr)
krt7_mask         = positive_any_gene_mask(expression_for_genes(adata, ["Krt7"]))
krt8_mask         = positive_any_gene_mask(expression_for_genes(adata, ["Krt8"]))
luminal_tf_mask   = positive_any_gene_mask(luminal_expr)
intermediate_mask = epcam_mask & krt7_mask & krt8_mask & luminal_tf_mask

kidney_low_mask  = score < KIDNEY_SCORE_MAX
strict_uro_mask  = (upk_mask | basal_mask | intermediate_mask) & kidney_low_mask

adata.obs["uro_umbrella_arm"]         = upk_mask
adata.obs["uro_basal_arm"]            = basal_mask
adata.obs["uro_intermediate_arm"]     = intermediate_mask
adata.obs["uro_low_kidney_epi_score"] = kidney_low_mask
adata.obs["strict_urothelium"]        = strict_uro_mask

print("\nUrothelial extraction gates:")
print_gate_count("any UPK > 0  (umbrella arm)", upk_mask, adata.n_obs)
print_gate_count(f"≥{BASAL_MIN_POSITIVE}/3 Krt5/Krt14/Trp63 > 0 (basal arm)", basal_mask, adata.n_obs)
print_gate_count("Epcam+Krt7+Krt8+Foxa1/Gata3 (intermediate arm)", intermediate_mask, adata.n_obs)
print_gate_count("any arm positive", upk_mask | basal_mask | intermediate_mask, adata.n_obs)
print_gate_count(f"{KIDNEY_SCORE_NAME} < {KIDNEY_SCORE_MAX}", kidney_low_mask, adata.n_obs)
print_gate_count("FINAL (any arm) AND low kidney score", strict_uro_mask, adata.n_obs)

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


###############################################################################
# Quick validation — arm breakdown in the extracted cells
###############################################################################

print("\nArm breakdown in extracted urothelial cells:")
print_gate_count("umbrella arm (any UPK)",
                 adata_uro.obs["uro_umbrella_arm"], adata_uro.n_obs)
print_gate_count(f"basal arm (≥{BASAL_MIN_POSITIVE}/3 Krt5/Krt14/Trp63)",
                 adata_uro.obs["uro_basal_arm"], adata_uro.n_obs)
print_gate_count("intermediate arm (Epcam + Krt7 +Krt8  Foxa1/Gata3)",
                 adata_uro.obs["uro_intermediate_arm"], adata_uro.n_obs)
