"""
04a_extract_uro_rawcounts_h5ad.py

Extract the true raw UMI count matrix for the 715 Uro-only cells (barcodes
written by 04_UrotheliumDevelopment_CytoTRACE.R's R companion) from the
source h5ad's `raw.X` slot. The RDS object used everywhere else in this
pipeline only carries a log-normalized `data` layer (see
[[project_chen2025_devatlas_structure]] memory note); this h5ad's `raw.X`
is confirmed to hold genuine integer counts (same shape/gene order as `X`),
so we pull the small Uro-only submatrix directly here rather than
re-preprocessing the whole 203k-cell object.

Reads raw/X directly via h5py (CSR format) instead of loading the full
anndata object, since we only need 715 of 203,139 rows.

Input:  UsedSingleCellsRawResults/RenalUrothelium/
        MultiOmicSpatialMouseKidney_SexSpecific_lifespan_Chen2025NatGenet.h5ad
        UrotheliumDevelopmentScripts/output/UrotheliumOnly_barcodes.txt
Output: UrotheliumDevelopmentScripts/output/UrotheliumOnly_rawcounts.mtx
        UrotheliumDevelopmentScripts/output/UrotheliumOnly_rawcounts_genes.txt
        UrotheliumDevelopmentScripts/output/UrotheliumOnly_rawcounts_barcodes.txt
"""

import h5py
import numpy as np
import scipy.sparse as sp
import scipy.io as sio

H5AD = ("/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/"
        "UsedSingleCellsRawResults/RenalUrothelium/"
        "MultiOmicSpatialMouseKidney_SexSpecific_lifespan_Chen2025NatGenet.h5ad")
OUT_DIR = ("/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/"
           "UsedSingleCells/UrotheliumDevelopmentScripts/output")

barcode_file = f"{OUT_DIR}/UrotheliumOnly_barcodes.txt"
with open(barcode_file) as fh:
    target_barcodes = [line.strip() for line in fh if line.strip()]
print(f"Loaded {len(target_barcodes)} target Uro barcodes")

with h5py.File(H5AD, "r") as f:
    all_barcodes = f["obs"]["_index"][:].astype(str)
    barcode_to_row = {bc: i for i, bc in enumerate(all_barcodes)}

    missing = [bc for bc in target_barcodes if bc not in barcode_to_row]
    if missing:
        raise SystemExit(f"{len(missing)} barcodes not found in h5ad obs, e.g. {missing[:5]}")

    row_idx = np.array([barcode_to_row[bc] for bc in target_barcodes], dtype=np.int64)

    indptr = f["raw"]["X"]["indptr"][:]
    n_genes = f["raw"]["var"]["_index"].shape[0]
    gene_ids = f["raw"]["var"]["_index"][:].astype(str)

    print(f"raw/X: {len(indptr) - 1} cells x {n_genes} genes total; extracting {len(row_idx)} rows ...")

    # CSR row-slicing: for each target row, pull its [start:end) slice of
    # data/indices directly rather than materializing the full 203k x 31.7k
    # matrix in memory.
    data_list, indices_list, new_indptr = [], [], [0]
    data_ds = f["raw"]["X"]["data"]
    indices_ds = f["raw"]["X"]["indices"]
    for r in row_idx:
        start, end = indptr[r], indptr[r + 1]
        data_list.append(data_ds[start:end])
        indices_list.append(indices_ds[start:end])
        new_indptr.append(new_indptr[-1] + (end - start))

    data = np.concatenate(data_list)
    indices = np.concatenate(indices_list)
    new_indptr = np.array(new_indptr, dtype=np.int64)

mat = sp.csr_matrix((data, indices, new_indptr), shape=(len(row_idx), n_genes))
print(f"Extracted matrix: {mat.shape}, nnz = {mat.nnz}, max value = {mat.max()}, "
      f"all-integer = {np.allclose(mat.data, np.round(mat.data))}")

# Write genes x cells (Seurat convention) as MatrixMarket + barcode/gene lists.
sio.mmwrite(f"{OUT_DIR}/UrotheliumOnly_rawcounts.mtx", mat.T.tocsr())
with open(f"{OUT_DIR}/UrotheliumOnly_rawcounts_genes.txt", "w") as fh:
    fh.write("\n".join(gene_ids) + "\n")
with open(f"{OUT_DIR}/UrotheliumOnly_rawcounts_barcodes.txt", "w") as fh:
    fh.write("\n".join(target_barcodes) + "\n")

print("Saved: UrotheliumOnly_rawcounts.mtx (+ genes.txt, barcodes.txt)")
