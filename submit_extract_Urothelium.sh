#!/bin/bash
#SBATCH --job-name=UUO_Urothelium
#SBATCH --output=/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/logs/UUO_Urothelium_%j.out
#SBATCH --error=/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/logs/UUO_Urothelium_%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=256G
#SBATCH --partition=himem

SCRIPT_DIR="/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"

mkdir -p "${SCRIPT_DIR}/logs"

# ---------------------------------------------------------------------------
# Step 1: Extract urothelial cells → Urothelium_cells.h5ad  (Python)
# ---------------------------------------------------------------------------
source /home/gdbecknelllab/xxw004/Software/miniconda3/etc/profile.d/conda.sh
conda activate cell2loc_env

echo "Python:      $(python --version)"
echo "Hostname:    $(hostname)"
echo "Start time:  $(date)"
echo ""

python "${SCRIPT_DIR}/extract_Urothelium_cells.py"

echo ""
echo "Python step done: $(date)"
echo ""

# ---------------------------------------------------------------------------
# Step 2: Convert h5ad → Seurat rds  (R)
# ---------------------------------------------------------------------------
conda deactivate

echo "Running R conversion ..."
Rscript "${SCRIPT_DIR}/convert_Urothelium_to_rds.R"

echo ""
echo "End time: $(date)"
