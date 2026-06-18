#!/usr/bin/env python3
"""
04b_scvi_AllUrothelium.py

scVI integration of all urothelium cells:
  bladder + kidney + ureter  ×  in-vivo + organoid  ×  multi-chemistry.

Input:  output/AllUrothelium_scvi_input/   (written by 04a_export_for_scvi.R)
Output: output/AllUrothelium_scvi/
  AllUrothelium_scvi.h5ad              – AnnData with latent + UMAP + clusters
  model/                               – saved scVI model (reloadable)
  AllUrothelium_scvi_clusters.csv      – cluster assignments (for R downstream)
  figures/                             – UMAP PDFs

Why these scVI settings
───────────────────────
• batch_key = "sample_id"
    Most granular batch variable.  scVI corrects for all sample-level effects
    (chemistry, lab, passage, etc.) without needing to know the hierarchy.

• categorical_covariate_keys = ["technology", "source"]
    Explicitly tells the encoder about two known confounders:
      technology — 10X / PIPseq / sci-RNA-seq3 chemistry differences
      source     — scVI-processed in-vivo vs raw-count organoid QC pipeline
    Harmony was unable to fully remove these because they dominate the top PCs;
    scVI conditions the encoder on them directly in a non-linear latent space.

• gene_likelihood = "normal"
    The input is log-normalised (R data layer), not raw UMI counts.
    scVI-source in-vivo files have scVI log-norm in their counts slot —
    original raw counts are unavailable.  The Gaussian (normal) decoder is
    the documented approach for pre-normalised data.

• n_latent = 30, n_layers = 2, max_epochs = 400
    Standard for a complex multi-source atlas of this size (~200k cells).
"""

import os
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

import pandas as pd
from scipy.io import mmread

import scanpy as sc
import scvi

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR  = "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/FinalUrotheliumCells"
SCR_DIR   = os.path.join(BASE_DIR, "UrotheliumScripts")
OUT_DIR   = os.path.join(SCR_DIR, "output")
SCVI_IN   = os.path.join(OUT_DIR, "AllUrothelium_scvi_input")
SCVI_OUT  = os.path.join(OUT_DIR, "AllUrothelium_scvi")
FIG_DIR   = os.path.join(SCVI_OUT, "figures")
MODEL_DIR = os.path.join(SCVI_OUT, "model")

for d in [SCVI_OUT, FIG_DIR, MODEL_DIR]:
    os.makedirs(d, exist_ok=True)

sc.settings.figdir = FIG_DIR

# ── Parameters ─────────────────────────────────────────────────────────────────
BATCH_KEY        = "sample_id"
CAT_COV_KEYS     = ["technology", "Categories"]
CONT_COV_KEYS    = ["pct_mt"]

N_LATENT         = 30
N_LAYERS         = 2
N_HIDDEN         = 128
GENE_LIKELIHOOD  = "normal"   # log-normalised input — must NOT be "nb" / "zinb"

MAX_EPOCHS       = 400
EARLY_STOP_PAT   = 20         # stop if ELBO doesn't improve for N epochs
LR               = 1e-3

N_NEIGHBORS      = 30
UMAP_MIN_DIST    = 0.3
LEIDEN_RESOLS    = [0.3, 0.5, 0.8, 1.0]
SCVI_SEED        = 0

COLOR_GROUPS = [
    "leiden_r0.5", "seurat_clusters",
    "tissue", "condition", "FinalConditionL1", "FinalConditionL2",
    "sample_id", "technology", "Categories", "paper", "source_GEO",
]


# ── Reproducibility ────────────────────────────────────────────────────────────
scvi.settings.seed = SCVI_SEED
sc.settings.verbosity = 2


# ── 1. Load exported data ──────────────────────────────────────────────────────
print("=" * 60)
print("Loading exported sparse matrix ...")
print("=" * 60)

mat      = mmread(os.path.join(SCVI_IN, "matrix.mtx")).tocsr().T   # cells × genes
barcodes = pd.read_csv(
    os.path.join(SCVI_IN, "barcodes.tsv"), header=None, names=["bc"]
)["bc"].tolist()
features = pd.read_csv(
    os.path.join(SCVI_IN, "features.tsv"), header=None, names=["gene"]
)["gene"].tolist()
meta = pd.read_csv(os.path.join(SCVI_IN, "metadata.csv"), index_col=0)

# Align metadata rows to barcode order
meta = meta.reindex(barcodes)

adata = sc.AnnData(
    X   = mat,
    obs = meta,
    var = pd.DataFrame(index=features),
)
adata.obs_names = barcodes
adata.var_names = features
adata.var_names_make_unique()

print(f"  {adata.n_obs:,} cells  ×  {adata.n_vars:,} genes")

# Store the log-normalised values in a named layer (keeps X clean for scanpy)
adata.layers["logdata"] = adata.X.copy()


# ── 2. Covariate hygiene ───────────────────────────────────────────────────────
print("\nChecking covariates ...")

for col in [BATCH_KEY] + CAT_COV_KEYS:
    if col not in adata.obs.columns:
        print(f"  WARNING: '{col}' not found — adding 'Unknown' placeholder")
        adata.obs[col] = "Unknown"
    n_na = adata.obs[col].isna().sum()
    if n_na > 0:
        print(f"  Filling {n_na} NA in '{col}' → 'Unknown'")
        adata.obs[col] = adata.obs[col].fillna("Unknown")
    adata.obs[col] = adata.obs[col].astype(str)
    print(f"  {col}: {adata.obs[col].nunique()} levels")

for col in CONT_COV_KEYS:
    if col not in adata.obs.columns:
        print(f"  WARNING: '{col}' not found — adding 0.0 placeholder")
        adata.obs[col] = 0.0
    n_na = adata.obs[col].isna().sum()
    if n_na > 0:
        med = adata.obs[col].median()
        adata.obs[col] = adata.obs[col].fillna(med)
    adata.obs[col] = adata.obs[col].astype(float)

# Drop covariates with only one level (scVI would error or ignore them anyway)
cat_covs_use  = [c for c in CAT_COV_KEYS  if adata.obs[c].nunique() > 1]
cont_covs_use = [c for c in CONT_COV_KEYS if c in adata.obs.columns]

dropped = set(CAT_COV_KEYS) - set(cat_covs_use)
if dropped:
    print(f"  Dropping single-level covariates: {dropped}")


# ── 3. scVI model setup ────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("Setting up scVI ...")
print("=" * 60)
print(f"  batch_key               : {BATCH_KEY}")
print(f"  categorical_covariates  : {cat_covs_use}")
print(f"  continuous_covariates   : {cont_covs_use}")
print(f"  gene_likelihood         : {GENE_LIKELIHOOD}")
print(f"  n_latent / n_layers     : {N_LATENT} / {N_LAYERS}")

scvi.model.SCVI.setup_anndata(
    adata,
    layer                      = "logdata",
    batch_key                  = BATCH_KEY,
    categorical_covariate_keys = cat_covs_use  if cat_covs_use  else None,
    continuous_covariate_keys  = cont_covs_use if cont_covs_use else None,
)

model = scvi.model.SCVI(
    adata,
    n_latent        = N_LATENT,
    n_layers        = N_LAYERS,
    n_hidden        = N_HIDDEN,
    gene_likelihood = GENE_LIKELIHOOD,
    dispersion      = "gene",
)
model.view_anndata_setup()
print(model)


# ── 4. Train ───────────────────────────────────────────────────────────────────
print(f"\nTraining (max {MAX_EPOCHS} epochs, early-stop patience {EARLY_STOP_PAT}) ...")
model.train(
    max_epochs              = MAX_EPOCHS,
    early_stopping          = True,
    early_stopping_patience = EARLY_STOP_PAT,
    plan_kwargs             = {"lr": LR},
)

# Plot training history
try:
    train_elbo = model.history["elbo_train"]
    val_elbo   = model.history.get("elbo_validation", None)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(train_elbo, label="train ELBO")
    if val_elbo is not None:
        ax.plot(val_elbo, label="val ELBO")
    ax.set_xlabel("Epoch"); ax.set_ylabel("ELBO"); ax.legend()
    ax.set_title("scVI training history")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "scvi_training_elbo.pdf"))
    plt.close(fig)
except Exception:
    pass

model.save(MODEL_DIR, overwrite=True)
print(f"  Model saved → {MODEL_DIR}")


# ── 5. Latent representation ───────────────────────────────────────────────────
print("\nExtracting latent representation ...")
adata.obsm["X_scVI"] = model.get_latent_representation()


# ── 6. kNN graph + UMAP ────────────────────────────────────────────────────────
print(f"Building kNN graph (k = {N_NEIGHBORS}, rep = X_scVI) ...")
sc.pp.neighbors(
    adata, use_rep="X_scVI",
    n_neighbors=N_NEIGHBORS, random_state=SCVI_SEED,
)

print("Running UMAP ...")
sc.tl.umap(adata, min_dist=UMAP_MIN_DIST, random_state=SCVI_SEED)


# ── 7. Leiden clustering ───────────────────────────────────────────────────────
print("Leiden clustering ...")
for res in LEIDEN_RESOLS:
    key = f"leiden_r{res}"
    sc.tl.leiden(adata, resolution=res, key_added=key, random_state=SCVI_SEED)
    n = adata.obs[key].nunique()
    print(f"  resolution {res:4.1f}  →  {n:2d} clusters  (obs['{key}'])")


# ── 8. UMAP figures ────────────────────────────────────────────────────────────
print("\nSaving UMAP figures ...")
for grp in COLOR_GROUPS:
    if grp not in adata.obs.columns:
        print(f"  Skipping '{grp}' (column not found)")
        continue
    n_lvl  = adata.obs[grp].nunique()
    legend = "on data" if n_lvl <= 25 else "right margin"
    try:
        sc.pl.umap(
            adata,
            color      = grp,
            legend_loc = legend,
            frameon    = False,
            title      = f"scVI UMAP  —  {grp}",
            save       = f"_scvi_{grp}.pdf",
            show       = False,
        )
        print(f"  Saved: umap_scvi_{grp}.pdf")
    except Exception as exc:
        print(f"  WARNING: plot for '{grp}' failed: {exc}")

# One combined overview: all Leiden resolutions on one page
fig, axes = plt.subplots(2, 2, figsize=(16, 14))
axes = axes.flatten()
for ax, res in zip(axes, LEIDEN_RESOLS):
    key = f"leiden_r{res}"
    sc.pl.umap(adata, color=key, ax=ax, frameon=False,
               title=f"Leiden r={res}", show=False, legend_loc="on data")
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "umap_scvi_leiden_overview.pdf"))
plt.close(fig)
print("  Saved: umap_scvi_leiden_overview.pdf")


# ── 9. Save outputs ────────────────────────────────────────────────────────────
h5ad_path = os.path.join(SCVI_OUT, "AllUrothelium_scvi.h5ad")
print(f"\nSaving h5ad → {h5ad_path}")
adata.write_h5ad(h5ad_path, compression="gzip")

# Cluster assignments CSV for downstream R analysis
clust_cols = [f"leiden_r{r}" for r in LEIDEN_RESOLS]
clust_cols = [c for c in clust_cols if c in adata.obs.columns]
adata.obs[clust_cols].to_csv(
    os.path.join(SCVI_OUT, "AllUrothelium_scvi_clusters.csv")
)
print("  Cluster assignments saved.")


# ── 10. Summary ────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("scVI integration complete")
print("=" * 60)
print(f"  Cells              : {adata.n_obs:,}")
print(f"  Genes              : {adata.n_vars:,}")
print(f"  Latent dimensions  : {N_LATENT}")
print(f"  Epochs trained     : {len(model.history['elbo_train'])}")
print(f"  Clusters (r=0.5)   : {adata.obs['leiden_r0.5'].nunique()}")
print(f"  Output dir         : {SCVI_OUT}")
