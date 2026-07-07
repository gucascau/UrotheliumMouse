#!/bin/bash
#SBATCH --job-name=test_plot_fixes
#SBATCH --output=/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/SpatialScripts/logs/test_plot_fixes_%j.out
#SBATCH --error=/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/SpatialScripts/logs/test_plot_fixes_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G
#SBATCH --time=00:20:00
#SBATCH --partition=himem

set -euo pipefail

SCRIPT_DIR="/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/SpatialScripts"
mkdir -p "${SCRIPT_DIR}/logs"

module purge
module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0

Rscript "${SCRIPT_DIR}/test_plot_fixes.R"

echo "Job finished: $(date)"
