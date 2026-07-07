#!/bin/bash
#SBATCH --job-name=visiumhd_uro_subsets
#SBATCH --output=/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/SpatialScripts/logs/visiumhd_uro_subsets_%j.out
#SBATCH --error=/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/SpatialScripts/logs/visiumhd_uro_subsets_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=04:00:00
#SBATCH --partition=himem

set -euo pipefail

SCRIPT_DIR="/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/SpatialScripts"
mkdir -p "${SCRIPT_DIR}/logs"
mkdir -p "${SCRIPT_DIR}/output"

module purge
module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0

Rscript "${SCRIPT_DIR}/04_VisiumHD_Urothelium_RegionSubsets.R"

echo "Job finished: $(date)"
