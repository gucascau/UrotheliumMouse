#!/usr/bin/env python3
"""
annotation_Lake.py — Annotate RenalUrothelium_preprocessed.h5ad using the
Lake 2025 snRNA-seq reference (mouse_kidney_snRNAseq_Lake2025_bioRxiv_V2.h5ad).

Pipeline:
  1. Load Lake reference; convert Ensembl IDs → gene symbols
  2. Extract raw counts from raw.X (CellXGene format)
  3. Load query (RenalUrothelium_preprocessed.h5ad); use counts layer
  4. Subset both to shared genes
  5. Train scVI on reference
  6. Concatenate ref + query; train scANVI (semi-supervised)
  7. Predict cell types on query cells; save CSV + annotated h5ad
  8. Diagnostic plots

Reference label: SubclassLevel1
  PT, TAL, PC, DCT, EC, IC, DTL, FIB, Myeloid, Lymphoid,
  CNT, PapE, POD, ATL, PEC, VSM/P, Ad, NEU

Author: Xin Wang
Date:   2026-04-24
"""

import os
import re
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scvi
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scipy.sparse as sp
import torch

# ── Reproducibility ───────────────────────────────────────────────────────────
scvi.settings.seed = 0
sc.settings.verbosity = 1

print("scvi-tools version:", scvi.__version__)
print("GPU available:     ", torch.cuda.is_available())

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR   = "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
SCRIPT_DIR = os.path.join(BASE_DIR, "RenalUrotheliumScripts")

REF_PATH   = os.path.join(BASE_DIR, "RenalUrothelium",
                          "mouse_kidney_snRNAseq_Lake2025_bioRxiv_V2.h5ad")
QUERY_PATH = os.path.join(SCRIPT_DIR, "output",
                          "RenalUrothelium_preprocessed.h5ad")
OUT_DIR    = os.path.join(SCRIPT_DIR, "output", "annotation_Lake")
os.makedirs(OUT_DIR, exist_ok=True)

LABEL_KEY     = "cell_type"
LAKE_LABEL    = "SubclassLevel1"   # coarse annotation in Lake 2025
REF_NAME      = "Lake"
BATCH_KEY     = "sample_id"        # batch variable already in both datasets


# ── Helper: Ensembl IDs → gene symbols (CellXGene format) ────────────────────
def ensembl_to_symbol(adata, name=""):
    if not str(adata.var_names[0]).startswith("ENSMUS"):
        print(f"  [{name}] var_names already look like symbols — skip.")
        return adata
    if "feature_name" not in adata.var.columns:
        raise ValueError("No 'feature_name' column found in adata.var.")
    symbols = adata.var["feature_name"].astype(str).values.copy()
    symbols = np.array([re.sub(r'_ENSMUSG\d+$', '', s) for s in symbols])
    bad = (symbols == "nan") | (symbols == "") | (symbols == "None")
    if bad.sum():
        print(f"  [{name}] WARNING: {bad.sum()} genes kept as Ensembl IDs.")
        symbols[bad] = adata.var_names[bad]
    adata.var_names = symbols
    adata.var_names_make_unique()
    print(f"  [{name}] Converted {adata.n_vars:,} var_names to gene symbols.")
    return adata


# ── Helper: extract integer counts from CellXGene h5ad ───────────────────────
def extract_raw_counts(adata, name=""):
    """
    CellXGene format: X = log-norm, raw.X = integer counts over full gene set.
    Returns an AnnData with X = raw integer counts (csr float32).
    """
    if adata.raw is not None:
        raw = adata.raw.to_adata()
        raw = ensembl_to_symbol(raw, name + "/raw")
        x   = sp.csr_matrix(raw.X).astype(np.float32)
        print(f"  [{name}] Using raw.X — {raw.n_obs:,} × {raw.n_vars:,}")
        out = ad.AnnData(X=x, obs=adata.obs.copy(), var=raw.var.copy())
        return out
    # Fallback: use X if already integer
    x    = adata.X
    vals = x.data if sp.issparse(x) else x.ravel()
    vals = vals[vals != 0][:5000]
    if len(vals) and np.all(vals == np.floor(vals)):
        print(f"  [{name}] X appears to be raw integer counts — using directly.")
        return adata
    # Last resort: counts layer
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"]
        print(f"  [{name}] Using 'counts' layer.")
        return adata
    raise ValueError(f"[{name}] Cannot find raw integer counts.")


###############################################################################
# STEP 1: Load Lake reference
###############################################################################

print(f"\n{'='*60}")
print(f"  Loading Lake 2025 reference")
print(f"{'='*60}")
adata_ref = sc.read_h5ad(REF_PATH)
print(f"  {adata_ref.n_obs:,} cells × {adata_ref.n_vars:,} genes")
print(f"  {LAKE_LABEL} values:\n",
      adata_ref.obs[LAKE_LABEL].value_counts().to_string())

adata_ref = extract_raw_counts(adata_ref, "Lake")
adata_ref = ensembl_to_symbol(adata_ref, "Lake")
adata_ref.obs[LABEL_KEY] = adata_ref.obs[LAKE_LABEL].astype(str)

n_labeled = adata_ref.obs[LABEL_KEY].notna().sum()
print(f"  Labeled cells: {n_labeled:,} / {adata_ref.n_obs:,}")
if n_labeled == 0:
    raise ValueError(f"Column '{LAKE_LABEL}' is entirely empty in reference.")


###############################################################################
# STEP 2: Load query
###############################################################################

print(f"\n{'='*60}")
print(f"  Loading query: RenalUrothelium_preprocessed.h5ad")
print(f"{'='*60}")
adata_query = sc.read_h5ad(QUERY_PATH)
print(f"  {adata_query.n_obs:,} cells × {adata_query.n_vars:,} genes")

# Use the raw integer counts layer from preprocessing
if "counts" in adata_query.layers:
    adata_query.X = adata_query.layers["counts"].copy()
    print("  Using 'counts' layer as X for query.")


###############################################################################
# STEP 3: Shared genes
###############################################################################

shared = adata_ref.var_names.intersection(adata_query.var_names)
print(f"\n  Shared genes (Lake ∩ query): {len(shared):,}")
if len(shared) < 200:
    raise ValueError(f"Only {len(shared)} shared genes — check gene name harmonisation.")

ref_sub   = adata_ref[:, shared].copy()
query_sub = adata_query[:, shared].copy()
print(f"  ref_sub  : {ref_sub.n_obs:,} × {ref_sub.n_vars:,}")
print(f"  query_sub: {query_sub.n_obs:,} × {query_sub.n_vars:,}")


###############################################################################
# STEP 4: scVI on reference
###############################################################################

print(f"\n{'='*60}")
print("  Training scVI on Lake reference")
print(f"{'='*60}")

scvi.model.SCVI.setup_anndata(ref_sub, labels_key=LABEL_KEY)
scvi_model = scvi.model.SCVI(ref_sub, n_latent=30)
scvi_model.train()
scvi_model.save(os.path.join(OUT_DIR, f"scvi_model_{REF_NAME}"), overwrite=True)
print("  scVI model saved.")

ref_sub.obsm["X_scVI"] = scvi_model.get_latent_representation()
sc.pp.neighbors(ref_sub, use_rep="X_scVI")
sc.tl.umap(ref_sub)
fig, ax = plt.subplots(figsize=(10, 7))
sc.pl.umap(ref_sub, color=LABEL_KEY,
           title=f"{REF_NAME} reference — {LAKE_LABEL}", ax=ax, show=False)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, f"{REF_NAME}_scVI_UMAP.pdf"), bbox_inches="tight")
plt.close()
print("  Saved: Lake_scVI_UMAP.pdf")


###############################################################################
# STEP 5: Prepare query + concatenate
###############################################################################

scvi.model.SCVI.prepare_query_anndata(query_sub, scvi_model)
query_sub.obs[LABEL_KEY] = pd.Categorical(["Unknown"] * query_sub.n_obs)

combined = ad.concat(
    [ref_sub, query_sub],
    label      = "batch",
    keys       = ["ref", "query"],
    index_unique = "-",
    merge      = "first",
)
print(f"\n  Combined: {combined.n_obs:,} cells × {combined.n_vars:,} genes")


###############################################################################
# STEP 6: scANVI
###############################################################################

print(f"\n{'='*60}")
print("  Training scANVI (semi-supervised)")
print(f"{'='*60}")

scvi.model.SCANVI.setup_anndata(
    combined,
    labels_key          = LABEL_KEY,
    batch_key           = "batch",
    unlabeled_category  = "Unknown",
)
scanvi_model = scvi.model.SCANVI(combined, n_latent=30)
scanvi_model.train(max_epochs=20)
scanvi_model.save(os.path.join(OUT_DIR, f"scanvi_model_{REF_NAME}"), overwrite=True)
print("  scANVI model saved.")


###############################################################################
# STEP 7: Predict on query
###############################################################################

query_mask  = combined.obs["batch"] == "query"
soft_preds  = scanvi_model.predict(combined[query_mask], soft=True)
label_names = scanvi_model.adata_manager.get_state_registry(
    scvi.REGISTRY_KEYS.LABELS_KEY)["categorical_mapping"]

prob_df     = pd.DataFrame(
    soft_preds,
    index   = combined.obs_names[query_mask],
    columns = label_names,
)
pred_labels = prob_df.idxmax(axis=1)
pred_conf   = prob_df.max(axis=1)

combined.obs.loc[query_mask, f"scanvi_{REF_NAME}_label"]      = pred_labels.values
combined.obs.loc[query_mask, f"scanvi_{REF_NAME}_confidence"] = pred_conf.values

# Confidence histogram
pred_conf.hist(bins=50, figsize=(7, 4))
plt.title(f"scANVI confidence — {REF_NAME} reference")
plt.xlabel("Confidence")
plt.ylabel("Cell count")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, f"scanvi_{REF_NAME}_confidence.pdf"))
plt.close()

# Combined UMAP
combined.obsm["X_scANVI"] = scanvi_model.get_latent_representation()
sc.pp.neighbors(combined, use_rep="X_scANVI")
sc.tl.umap(combined)
fig, axes = plt.subplots(1, 3, figsize=(24, 7))
for ax, col, title in zip(axes,
                           [LABEL_KEY, "batch", f"scanvi_{REF_NAME}_label"],
                           ["Cell type (ref label)",
                            "Reference vs query",
                            f"scANVI predicted ({REF_NAME})"]):
    sc.pl.umap(combined, color=col, title=title, ax=ax, show=False)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, f"{REF_NAME}_combined_scANVI_UMAP.pdf"),
            bbox_inches="tight")
plt.close()
print("  Saved: Lake_combined_scANVI_UMAP.pdf")

# Label distribution
print(f"\n  Predicted label distribution (query):")
print(pred_labels.value_counts().to_string())
print(f"\n  Mean confidence: {pred_conf.mean():.3f}")
print(f"  Cells ≥ 0.9 confidence: "
      f"{(pred_conf >= 0.9).sum():,} / {query_mask.sum():,}")


###############################################################################
# STEP 8: Save predictions CSV + annotated query h5ad
###############################################################################

extra_cols = [c for c in adata_query.obs.columns
              if c in ("condition", "source", "sample_id", "technology",
                       "cell_type_original", "epithelial_compartment",
                       "assign_method")]

preds = combined.obs.loc[
    query_mask,
    [f"scanvi_{REF_NAME}_label", f"scanvi_{REF_NAME}_confidence"],
].copy()
if extra_cols:
    preds = preds.join(adata_query.obs[extra_cols], how="left")

csv_path = os.path.join(OUT_DIR, f"RenalUrothelium_scanvi_{REF_NAME}_annotations.csv")
preds.to_csv(csv_path)
print(f"\n  Predictions saved → {csv_path}")

# Annotate and save the query h5ad
adata_query.obs[f"scanvi_{REF_NAME}_label"]      = pred_labels.values
adata_query.obs[f"scanvi_{REF_NAME}_confidence"] = pred_conf.values
h5ad_path = os.path.join(OUT_DIR,
                         f"RenalUrothelium_annotated_{REF_NAME}.h5ad")
adata_query.write_h5ad(h5ad_path)
print(f"  Annotated query saved → {h5ad_path}")

print("\n===== Lake annotation complete =====")
