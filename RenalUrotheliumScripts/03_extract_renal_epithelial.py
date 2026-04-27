#!/usr/bin/env python3
"""
03_extract_renal_epithelial.py — Subset the preprocessed AnnData to renal
epithelial cells using a two-pass hybrid strategy:

  Pass 1 (annotation)  — match cell_type_original against known label sets
  Pass 2 (markers)     — for cells without a recognised label, score against
                         per-compartment marker gene panels and assign any
                         cell whose top-scoring compartment exceeds MIN_SCORE

Four compartments:
  Tubule               — proximal tubule, loop of Henle, distal convoluted
                         tubule, connecting tubule, nephron progenitors
  Urothelium           — urothelial / transitional epithelium
  Collecting_duct      — principal cells, intercalated cells A/B
  Glomerular_epithelial— podocytes, parietal epithelial cells (PEC)

Outputs
-------
  output/RenalUrothelium_renal_epithelial.h5ad
    obs column  epithelial_compartment : compartment name
    obs column  assign_method          : "annotation" | "marker_score"
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR   = "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
SCRIPT_DIR = os.path.join(BASE_DIR, "RenalUrotheliumScripts")
OUT_DIR    = os.path.join(SCRIPT_DIR, "output")

IN_PATH  = os.path.join(OUT_DIR, "RenalUrothelium_integrated.h5ad")
OUT_PATH = os.path.join(OUT_DIR, "RenalUrothelium_renal_epithelial.h5ad")

# ── Score threshold for marker-based pass ────────────────────────────────────
# sc.tl.score_genes returns a z-score–like value; 0.0 is a conservative floor
MIN_SCORE = 0.0

# ── Annotation labels per compartment ────────────────────────────────────────
COMPARTMENT_LABELS = {
    "Tubule": {
        # Proximal tubule — all segment nomenclatures across datasets
        "PTS1", "PTS2", "PTS3", "PTS3T2",
        "PT(S1)", "PT(S2)", "PT(S3)", "PT(S1-S2)", "PT",
        "Dediff. PT", "Prolif. PT",
        # Distal convoluted tubule
        "DCT", "DCT-CNT",
        # Connecting tubule
        "CNT",
        # Loop of Henle (all segments)
        "LOH", "LOH_AL", "LOH_DL", "LOH_AL_proliferating",
        "MTAL", "CTAL", "TAL",
        "LH(AL)", "LH(DL)",
        "DTL", "DTL-ATL", "ATL",
        # Distal tubule shorthand used in some datasets
        "D1", "D2",
        # Macula densa (specialised TAL epithelial)
        "MD",
        # Nephron / ureteric-bud progenitors
        "NP", "NP_proliferate", "UBP",
    },
    "Urothelium": {
        "Urotherlium",   # typo present in source dataset
        "Uro",
    },
    "Collecting_duct": {
        "CD_PC", "CD-PC", "CD_PC_Mix", "CD-Trans",
        "CD_IC",
        "IC", "IC-A", "IC-B", "ICA", "ICB",
        "PC",
    },
    "Glomerular_epithelial": {
        "Podo", "Pod",   # podocytes
        "PEC",           # parietal epithelial cells
    },
}

# ── Marker gene panels per compartment (mouse gene symbols) ──────────────────
COMPARTMENT_MARKERS = {
    "Tubule": [
        # Proximal tubule
        "Lrp2", "Slc34a1", "Slc13a3", "Cubn", "Hnf4a",
        "Slc27a2", "Slc22a6", "Slc22a8", "Fxyd2",
        # Loop of Henle / TAL
        "Umod", "Slc12a1", "Cldn16", "Cldn10",
        # Distal convoluted tubule
        "Slc12a3", "Pvalb", "Egf",
        # Connecting tubule / nephron progenitor
        "Calb1", "Cited1", "Six2",
    ],
    "Urothelium": [
        "Upk1a", "Upk1b", "Upk2", "Upk3a",
        "Krt5", "Krt14", "Krt20", "Krt8", "Krt18",
        "Trp63", "Upk3b",
    ],
    "Collecting_duct": [
        # Principal cells
        "Aqp2", "Avpr2", "Hsd11b2", "Scnn1a",
        # Intercalated cells
        "Atp6v1b1", "Atp6v0d2", "Slc4a1", "Foxi1",
        "Slc26a4", "Hmx2",
    ],
    "Glomerular_epithelial": [
        # Podocytes
        "Nphs1", "Nphs2", "Podxl", "Wt1", "Synpo",
        # Parietal epithelial cells
        "Pax8", "Cldn1", "Akr1b7", "Ptpro",
    ],
}

COMPARTMENT_ORDER = ["Tubule", "Urothelium",
                     "Collecting_duct", "Glomerular_epithelial"]

# Reverse lookup: label → compartment
LABEL_TO_COMPARTMENT = {
    lbl: comp
    for comp, labels in COMPARTMENT_LABELS.items()
    for lbl in labels
}


###############################################################################
# Load
###############################################################################

print(f"Loading {IN_PATH} ...")
adata = sc.read_h5ad(IN_PATH)
print(f"  Total: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

if "cell_type_original" not in adata.obs.columns:
    raise KeyError("'cell_type_original' column not found in adata.obs")

ct = adata.obs["cell_type_original"].astype(str)


###############################################################################
# Pass 1 — annotation-based assignment
###############################################################################

print("\n=== Pass 1: annotation-based ===")
compartment_col = ct.map(LABEL_TO_COMPARTMENT)   # NaN for unmatched cells
n_annot = compartment_col.notna().sum()
print(f"  Cells assigned by annotation: {n_annot:,}")


###############################################################################
# Pass 2 — marker-based assignment for unannotated cells
###############################################################################

print("\n=== Pass 2: marker-based (unannotated cells) ===")
unannotated_mask = compartment_col.isna()
n_unanno = unannotated_mask.sum()
print(f"  Unannotated cells to score  : {n_unanno:,}")

# Subset to unannotated cells for scoring
adata_un = adata[unannotated_mask].copy()

score_df = pd.DataFrame(index=adata_un.obs_names)

for comp in COMPARTMENT_ORDER:
    genes = [g for g in COMPARTMENT_MARKERS[comp] if g in adata_un.var_names]
    missing = [g for g in COMPARTMENT_MARKERS[comp] if g not in adata_un.var_names]
    print(f"  {comp}: {len(genes)} / {len(COMPARTMENT_MARKERS[comp])} markers "
          f"present  (absent: {missing if missing else 'none'})")
    if genes:
        sc.tl.score_genes(adata_un, gene_list=genes, score_name=f"_score_{comp}")
        score_df[comp] = adata_un.obs[f"_score_{comp}"].values
    else:
        score_df[comp] = -np.inf

# Assign compartment to cells where the best score exceeds MIN_SCORE
best_score = score_df.max(axis=1)
best_comp  = score_df.idxmax(axis=1)

marker_compartment = pd.Series(index=adata_un.obs_names, dtype=object)
marker_compartment[best_score >= MIN_SCORE] = best_comp[best_score >= MIN_SCORE]

n_marker = marker_compartment.notna().sum()
print(f"\n  Cells newly assigned by markers: {n_marker:,}")
print("  Marker-assigned breakdown:")
print(marker_compartment.value_counts().to_string())


###############################################################################
# Merge both passes
###############################################################################

# Start with annotation pass; fill gaps with marker pass
merged = compartment_col.copy()
merged[unannotated_mask] = marker_compartment.values

# Record assignment method
method = pd.Series("", index=adata.obs_names)
method[compartment_col.notna()]                                    = "annotation"
method[unannotated_mask & pd.notna(pd.Series(marker_compartment.values,
                                              index=adata.obs_names[unannotated_mask]))] = "marker_score"

adata.obs["epithelial_compartment"] = merged
adata.obs["assign_method"]          = method

final_mask = adata.obs["epithelial_compartment"].notna()
n_total    = final_mask.sum()
print(f"\n  Total renal epithelial cells: {n_total:,} / {adata.n_obs:,} "
      f"({100 * n_total / adata.n_obs:.2f}%)")


###############################################################################
# Composition report
###############################################################################

print("\n  Cells per compartment (both passes):")
print(adata.obs.loc[final_mask, "epithelial_compartment"].value_counts().to_string())

print("\n  Assignment method breakdown:")
print(adata.obs.loc[final_mask, "assign_method"].value_counts().to_string())

print("\n  Cells per cell_type_original (selected):")
print(adata.obs.loc[final_mask, "cell_type_original"].value_counts().to_string())

if "condition" in adata.obs.columns:
    print("\n  Cells per condition:")
    print(adata.obs.loc[final_mask, "condition"].value_counts().to_string())

if "source" in adata.obs.columns:
    print("\n  Cells per source:")
    print(adata.obs.loc[final_mask, "source"].value_counts().to_string())


###############################################################################
# Subset and save
###############################################################################

adata_epi = adata[final_mask].copy()
print(f"\nSaving {adata_epi.n_obs:,} renal epithelial cells → {OUT_PATH}")
adata_epi.write_h5ad(OUT_PATH)
print("Done.")
print(f"\n===== Extraction complete =====")
print(f"  Cells : {adata_epi.n_obs:,}")
print(f"  Genes : {adata_epi.n_vars:,}")
print(f"  Output: {OUT_PATH}")
