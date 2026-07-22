#!/bin/bash
#SBATCH --job-name=visiumlow_spp1itgav1
#SBATCH --output=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/VisiumLowDevelopmentScripts/logs/visiumlow_spp1itgav1_%j.out
#SBATCH --error=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/VisiumLowDevelopmentScripts/logs/visiumlow_spp1itgav1_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --partition=himem

set -euo pipefail

SCRIPT_DIR="/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/VisiumLowDevelopmentScripts"
mkdir -p "${SCRIPT_DIR}/logs"
mkdir -p "${SCRIPT_DIR}/output"

module purge
module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0

echo "Hostname  : $(hostname)"
echo "Start time: $(date)"

Rscript "${SCRIPT_DIR}/07_VisiumLow_Spp1Itgav1_SpatialAnalysis.R"

echo "Job finished: $(date)"
