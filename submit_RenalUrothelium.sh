#!/bin/bash
#SBATCH --job-name=RenalUrothelium
#SBATCH --output=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/logs/RenalUrothelium_%j.out
#SBATCH --error=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/logs/RenalUrothelium_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=128G
#SBATCH --time=4:00:00
#SBATCH --partition=himem

SCRIPT_DIR="/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"

mkdir -p "${SCRIPT_DIR}/logs"

echo "Hostname:   $(hostname)"
echo "Start time: $(date)"
echo ""

# ── Step 1: Extract renal urothelial cells from scvi_qc_integrated.h5ad ───────
source /home/gdbecknelllab/xxw004/Software/miniconda3/etc/profile.d/conda.sh
conda activate cell2loc_env

echo "=== Step 1: Python extraction ==="
python "${SCRIPT_DIR}/extract_RenalUrothelium.py"
echo "Python done: $(date)"
echo ""

conda deactivate

# ── Step 2: Convert h5ad → Seurat rds ─────────────────────────────────────────
module purge
module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0

echo "=== Step 2: R conversion (h5ad → rds) ==="
Rscript "${SCRIPT_DIR}/convert_RenalUrothelium_to_rds.R"
echo "Conversion done: $(date)"
echo ""

# ── Step 3: Harmony reclustering ──────────────────────────────────────────────
echo "=== Step 3: Harmony reclustering ==="
Rscript "${SCRIPT_DIR}/cluster_RenalUrothelium_harmony.R"
echo "Reclustering done: $(date)"

echo ""
echo "End time: $(date)"

sacct -j "${SLURM_JOB_ID}" \
  --format=JobID,User,Start,End,Elapsed,TotalCPU,State,NodeList \
  --noheader
