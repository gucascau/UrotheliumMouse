#!/bin/bash
#SBATCH --job-name=RenalUro_fullgene
#SBATCH --output=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/RenalUrotheliumScripts/logs/fullgene_%j.out
#SBATCH --error=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/RenalUrotheliumScripts/logs/fullgene_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=1300G
#SBATCH --gres=gpu:1
#SBATCH --time=72:00:00
#SBATCH --partition=himem

# ── Full-gene pipeline ────────────────────────────────────────────────────────
#
#   Step 10 : Build full-gene h5ad from original per-sample files
#             (~25k genes, all 1.1M cells incl. reference atlases)
#   Step 11 : Attach scVI/scANVI embeddings from the existing 3,000-HVG models,
#             then store all genes + UMAP + clusters
#   Step 10b: Convert full-gene scVI-integrated h5ad → Seurat RDS
#
# Memory note: UUOProjectObject_F.h5ad is 44 GB on disk (~150 GB in RAM).
# Total peak RAM ~400–600 GB. 1 300 G requested for safety.
# GPU requested for faster scVI/scANVI latent extraction; Step 11 does not retrain.

SCRIPT_DIR="/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/RenalUrotheliumScripts"
OUT_DIR="${SCRIPT_DIR}/output"

mkdir -p "${SCRIPT_DIR}/logs" "${OUT_DIR}"

echo "======================================================"
echo "  RenalUrothelium Full-Gene Pipeline"
echo "  Host  : $(hostname)"
echo "  Start : $(date)"
echo "======================================================"

source /home/gdbecknelllab/xxw004/Software/miniconda3/etc/profile.d/conda.sh
conda activate cell2loc_env
export SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-16}

# ── Step 10: Build full-gene h5ad ─────────────────────────────────────────────
echo ""
echo "=== Step 10: Build full-gene h5ad (all cells, ~25k genes) ==="
echo "Start: $(date)"
python "${SCRIPT_DIR}/10_build_fullgene_h5ad.py"
if [ $? -ne 0 ]; then echo "ERROR in Step 10"; exit 1; fi
echo "Done:  $(date)"

# ── Step 11: scVI integration ─────────────────────────────────────────────────
echo ""
echo "=== Step 11: scVI/scANVI embeddings (existing 3000-HVG models, full gene object) ==="
echo "Start: $(date)"
python "${SCRIPT_DIR}/11_scvi_integrate_fullgene.py"
if [ $? -ne 0 ]; then echo "ERROR in Step 11"; exit 1; fi
echo "Done:  $(date)"

conda deactivate

# ── Step 10b: Convert to Seurat RDS ──────────────────────────────────────────
echo ""
echo "=== Step 10b: Convert full-gene scVI h5ad → Seurat RDS ==="
echo "Start: $(date)"

module purge
module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0
export LD_LIBRARY_PATH=/gpfs0/export/apps/easybuild/software/GCCcore/12.3.0/lib64:${LD_LIBRARY_PATH}

# Point 10b at the scVI-integrated output instead of the plain fullgene h5ad
# by temporarily symlinking — or just pass the path via env var
export FULLGENE_H5AD="${OUT_DIR}/RenalUrothelium_allcells_scvi.h5ad"
export FULLGENE_RDS="${OUT_DIR}/RenalUrothelium_allcells_scvi.rds"

Rscript - <<'REOF'
suppressPackageStartupMessages({
  library(zellkonverter); library(Seurat); library(SingleCellExperiment)
})
h5ad_in  <- Sys.getenv("FULLGENE_H5AD")
int_rds  <- file.path(dirname(h5ad_in), "RenalUrothelium_integrated_annotated.rds")
rds_out  <- Sys.getenv("FULLGENE_RDS")

cat("Reading", h5ad_in, "\n")
sce <- readH5AD(h5ad_in, use_hdf5 = FALSE)
cat("  SCE dims:", nrow(sce), "genes x", ncol(sce), "cells\n")

cat("Converting SCE -> Seurat ...\n")
so <- as.Seurat(sce, counts = "counts", data = "lognorm")
rm(sce); gc()

cat("Transferring reductions from integrated RDS ...\n")
so_int    <- readRDS(int_rds)
shared    <- intersect(colnames(so), colnames(so_int))
so_int    <- so_int[, shared]
so        <- so[,    shared]

for (red in c("scVI", "scANVI", "umap_scVI", "umap_scANVI")) {
  if (red %in% names(so_int@reductions)) {
    so[[paste0("orig_", red)]] <- CreateDimReducObject(
      embeddings = Embeddings(so_int, red),
      key        = paste0("orig", gsub("_","",red), "_"),
      assay      = DefaultAssay(so)
    )
    cat("  Transferred (original 3k-HVG):", red, "\n")
  }
}
rm(so_int); gc()

cat("Saving ->", rds_out, "\n")
saveRDS(so, rds_out)
cat("Done.\n  Cells:", ncol(so), "\n  Genes:", nrow(so), "\n")
REOF

if [ $? -ne 0 ]; then echo "ERROR in Step 10b"; exit 1; fi
echo "Done:  $(date)"

echo ""
echo "======================================================"
echo "  Full-gene pipeline complete: $(date)"
echo "======================================================"

sacct -j "${SLURM_JOB_ID}" \
  --format=JobID,User,Start,End,Elapsed,TotalCPU,State,NodeList \
  --noheader
