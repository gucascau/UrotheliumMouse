#!/bin/bash
#SBATCH --job-name=visiumhd_deconv
#SBATCH --output=/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/SpatialScripts/logs/visiumhd_deconv_%j.out
#SBATCH --error=/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/SpatialScripts/logs/visiumhd_deconv_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=800G
#SBATCH --time=48:00:00
#SBATCH --partition=himem

set -euo pipefail

SCRIPT_DIR="/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/SpatialScripts"
mkdir -p "${SCRIPT_DIR}/logs"
mkdir -p "${SCRIPT_DIR}/output"

module purge
module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0

# Rscript "${SCRIPT_DIR}/02_VisiumHD_Deconvolution.r"
Rscript "${SCRIPT_DIR}/02b_VisiumHD_Deconvolution_BothKidneys.R"

echo "Job finished: $(date)"
