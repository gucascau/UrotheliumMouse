################################################################################
# 00_export_Chen2025_rawcounts.py
#
# The converted Chen2025 Seurat RDS (MultiOmicSpatialMouseKidney_..._Chen2025
# NatGenet_zellkonvertedConverted.rds) only carries a log-normalized "data"
# layer -- no raw counts survived the zellkonverter h5ad->Seurat conversion.
# RCTD (spacexr) needs raw UMI counts + nUMI for its Poisson-based model, not
# normalized floats.
#
# The original h5ad (same directory, pre-conversion) still has raw counts in
# its CellxGene-schema `raw.X` slot (confirmed integer-valued, same 203,139
# cells x 31,671 genes, same cell/gene order as `X` and as the converted RDS
# -- checked directly against the RDS's colnames/rownames). This script
# exports raw.X to a standard 10x-style mtx/barcodes/features trio so the R
# deconvolution script (05_VisiumLow_Deconvolution.R) can load it with
# Seurat::Read10X() and attach it as counts, matched by barcode.
#
# One-time preprocessing step -- not part of the 01-04 numbered R pipeline,
# hence "00" and .py rather than .R (h5py is more direct than round-tripping
# through anndata2ri/reticulate for a single sparse-matrix export).
#
# Input:  RawMouseSingleCellDatasets/.../Chen2025NatGenet.h5ad
# Output: VisiumLowDevelopmentScripts/output/Chen2025_rawcounts/
#         {matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz}
################################################################################

import gzip
import h5py
import numpy as np
from scipy.io import mmwrite
from scipy.sparse import csr_matrix
import os

H5AD = "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCellsRawResults/RenalUrothelium/MultiOmicSpatialMouseKidney_SexSpecific_lifespan_Chen2025NatGenet.h5ad"
OUT_DIR = "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/VisiumLowDevelopmentScripts/output/Chen2025_rawcounts"
os.makedirs(OUT_DIR, exist_ok=True)

print("==> Reading raw/X (CSR, cells x genes) ...")
with h5py.File(H5AD, "r") as f:
    shape = tuple(f["raw"]["X"].attrs["shape"])
    data = f["raw"]["X"]["data"][:]
    indices = f["raw"]["X"]["indices"][:]
    indptr = f["raw"]["X"]["indptr"][:]
    barcodes = [b.decode() for b in f["obs"]["_index"][:]]
    gene_ids = [b.decode() for b in f["var"]["_index"][:]]
    gene_symbols = [b.decode() for b in f["var"]["gene_symbols"][:]]

print(f"  shape (cells x genes): {shape}, nnz: {len(data)}")
raw_csr = csr_matrix((data, indices, indptr), shape=shape)

# Sanity check: raw counts should be non-negative integers.
sample = raw_csr.data[:1000]
assert np.all(sample >= 0) and np.allclose(sample, np.round(sample)), \
    "raw/X does not look like integer counts -- aborting."

print("==> Transposing to genes x cells (10x/Read10X convention) ...")
counts_gxc = raw_csr.transpose().tocsr().astype(np.int32)

print("==> Writing matrix.mtx.gz ...")
mtx_path = os.path.join(OUT_DIR, "matrix.mtx")
mmwrite(mtx_path, counts_gxc, field="integer")
with open(mtx_path, "rb") as f_in, gzip.open(mtx_path + ".gz", "wb") as f_out:
    f_out.writelines(f_in)
os.remove(mtx_path)

print("==> Writing barcodes.tsv.gz and features.tsv.gz ...")
with gzip.open(os.path.join(OUT_DIR, "barcodes.tsv.gz"), "wt") as f:
    f.write("\n".join(barcodes) + "\n")

with gzip.open(os.path.join(OUT_DIR, "features.tsv.gz"), "wt") as f:
    for gid, gsym in zip(gene_ids, gene_symbols):
        f.write(f"{gid}\t{gsym}\tGene Expression\n")

print("==> Done. Wrote:", OUT_DIR)
