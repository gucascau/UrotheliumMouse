#!/usr/bin/env python3
"""
03_extract_renal_epithelial.py — Subset the preprocessed AnnData to renal
epithelial cells using a two-pass hybrid strategy:

  Pass 1 (annotation)  — use scANVI-predicted labels from two reference atlases
                         (scanvi_Lake_label, scanvi_MKA_label); a cell is
                         assigned if EITHER predicts an epithelial type
  Pass 2 (markers)     — for cells not predicted as epithelial by either
                         reference, score against per-compartment marker gene
                         panels and assign any cell whose top-scoring
                         compartment exceeds MIN_SCORE

Four compartments:
  Tubule               — proximal tubule, loop of Henle, distal convoluted
                         tubule, connecting tubule
  Urothelium           — urothelial / papillary epithelium
  Collecting_duct      — principal cells, intercalated cells A/B
  Glomerular_epithelial— podocytes, parietal epithelial cells (PEC)

Outputs
-------
  output/RenalUrothelium_renal_epithelial.h5ad
    obs column  epithelial_compartment : compartment name
    obs column  assign_method          : "annotation_Lake" | "annotation_MKA"
                                         | "annotation_both" | "marker_score"
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR   = "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
SCRIPT_DIR = os.path.join(BASE_DIR, "RenalUrotheliumScripts")
OUT_DIR    = os.path.join(SCRIPT_DIR, "output")

IN_PATH  = os.path.join(OUT_DIR, "RenalUrothelium_integrated_annotated.h5ad")
OUT_PATH = os.path.join(OUT_DIR, "RenalUrothelium_renal_epithelial.h5ad")

# ── Score threshold for marker-based pass ────────────────────────────────────
# sc.tl.score_genes returns a z-score–like value; 0.0 is a conservative floor
MIN_SCORE = 0.0

# ── Annotation label sets — curated from scanvi_Lake_label / scanvi_MKA_label ─
LAKE_COMPARTMENT_LABELS = {
    "Tubule": {
        "PT", "DCT", "DTL", "TAL", "ATL", "CNT", "PapE",
    },
    "Collecting_duct": {
        "PC", "IC",
    },
    "Glomerular_epithelial": {
        "PEC", "POD",
    },
}

MKA_COMPARTMENT_LABELS = {
    "Tubule": {
        "PTS1", "PTS2", "PTS3", "PTS3T2",
        "DCT", "DCT-CNT", "CNT",
        "MTAL", "CTAL",
        "DTL", "ATL", "DTL-ATL",
        "LOH", "MD",
    },
    "Collecting_duct": {
        "PC", "ICA", "ICB", "CD-Trans",
    },
    "Glomerular_epithelial": {
        "Podo", "PEC",
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

COMPARTMENT_ORDER = ["Tubule", "Urothelium", "Collecting_duct", "Glomerular_epithelial"]

LAKE_LABEL_TO_COMPARTMENT = {
    lbl: comp for comp, labels in LAKE_COMPARTMENT_LABELS.items() for lbl in labels
}
MKA_LABEL_TO_COMPARTMENT = {
    lbl: comp for comp, labels in MKA_COMPARTMENT_LABELS.items() for lbl in labels
}


###############################################################################
# Load
###############################################################################

print(f"Loading {IN_PATH} ...")
adata = sc.read_h5ad(IN_PATH)
print(f"  Total: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

for col in ("scanvi_Lake_label", "scanvi_MKA_label"):
    if col not in adata.obs.columns:
        raise KeyError(f"'{col}' column not found in adata.obs")

ct_lake = adata.obs["scanvi_Lake_label"].astype(str)
ct_mka  = adata.obs["scanvi_MKA_label"].astype(str)


###############################################################################
# Pass 1 — annotation-based assignment
###############################################################################

print("\n=== Pass 1: annotation-based (Lake + MKA scANVI predictions) ===")

lake_comp = ct_lake.map(LAKE_LABEL_TO_COMPARTMENT)   # NaN if not epithelial in Lake
mka_comp  = ct_mka.map(MKA_LABEL_TO_COMPARTMENT)    # NaN if not epithelial in MKA

# Include cell if EITHER annotation predicts epithelial; Lake takes priority on conflict
compartment_col = lake_comp.combine_first(mka_comp)

# Track which reference(s) drove the assignment
method_annot = pd.Series("", index=adata.obs_names)
method_annot[lake_comp.notna() & mka_comp.isna()]  = "annotation_Lake"
method_annot[lake_comp.isna()  & mka_comp.notna()] = "annotation_MKA"
method_annot[lake_comp.notna() & mka_comp.notna()]  = "annotation_both"

n_lake  = (method_annot == "annotation_Lake").sum()
n_mka   = (method_annot == "annotation_MKA").sum()
n_both  = (method_annot == "annotation_both").sum()
n_annot = compartment_col.notna().sum()
print(f"  Lake only        : {n_lake:,}")
print(f"  MKA only         : {n_mka:,}")
print(f"  Both (Lake wins) : {n_both:,}")
print(f"  Total assigned   : {n_annot:,}")


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

# Fill annotation gaps with marker-based assignments
merged = compartment_col.copy()
merged[unannotated_mask] = marker_compartment.values

# Final method column: annotation source or marker_score
method = method_annot.copy()
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

print("\n  scanvi_Lake_label breakdown (selected):")
print(adata.obs.loc[final_mask, "scanvi_Lake_label"].value_counts().to_string())

print("\n  scanvi_MKA_label breakdown (selected):")
print(adata.obs.loc[final_mask, "scanvi_MKA_label"].value_counts().to_string())

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
