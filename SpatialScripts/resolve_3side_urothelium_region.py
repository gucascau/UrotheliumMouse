################################################################################
# resolve_3side_urothelium_region.py
#
# Selected_3SideUrotheliumRegion_V2.csv contains VisiumHD *segmented-cell*
# barcodes ("cellid_XXXXXXXX-1"), exported from the Loupe Browser cell-
# segmentation view. These IDs do not exist in the Spatial.008um assay used
# by the integrated Seurat object (square-bin assay), so they must first be
# resolved to their constituent square_008um bin barcodes via Space Ranger's
# barcode_mappings.parquet.
#
# Verified sample of origin: kidney3p (Visium_HD_3prime_Mouse_Kidney).
# All 10,187 cell IDs in the CSV resolve against kidney3p's mapping table
# (0 unmatched); kidney5p's mapping table only resolves ~94% (a false-
# positive rate expected from coincidental ID overlap between samples).
#
# Requires pyarrow (not available in the project's default R/conda envs);
# run once from a throwaway venv, output is a plain-text bin list consumed
# by 04_VisiumHD_Urothelium_RegionSubsets.R — no parquet parsing needed in R.
################################################################################

import pyarrow.parquet as pq
import pyarrow.compute as pc

HD_BASE = "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSpatialData/VisiumHD"
CSV_IN  = f"{HD_BASE}/Selected_3SideUrotheliumRegion_V2.csv"
PARQUET = f"{HD_BASE}/Visium_HD_3prime_Mouse_Kidney/Visium_HD_3prime_Mouse_Kidney_barcode_mappings.parquet"
OUT     = "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/SpatialScripts/output/Selected_3SideUrotheliumRegion_V2_resolved_008um_bins.txt"

csv_ids = set()
with open(CSV_IN) as f:
    next(f)
    for line in f:
        csv_ids.add(line.strip().split(",")[0])

tbl = pq.read_table(PARQUET, columns=["square_008um", "cell_id", "in_cell"])
mask = pc.and_(pc.is_valid(tbl["cell_id"]), tbl["in_cell"])
sub = tbl.filter(mask)

cell_ids = sub["cell_id"].to_pylist()
bins008 = sub["square_008um"].to_pylist()

matched = len(set(c for c in cell_ids if c in csv_ids))
print(f"Matched unique cell_ids: {matched} / {len(csv_ids)}")

bin_set = sorted(set(b for c, b in zip(cell_ids, bins008) if c in csv_ids))
print(f"Resolved to {len(bin_set)} unique square_008um bins")

with open(OUT, "w") as out:
    for b in bin_set:
        out.write(b + "\n")
print("Saved:", OUT)
