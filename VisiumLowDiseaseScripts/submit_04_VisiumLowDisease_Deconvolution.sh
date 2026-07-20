#!/bin/bash
#SBATCH --job-name=visiumlowdisease_deconv
#SBATCH --output=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/VisiumLowDiseaseScripts/logs/visiumlowdisease_deconv_%j.out
#SBATCH --error=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/VisiumLowDiseaseScripts/logs/visiumlowdisease_deconv_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=1000G
#SBATCH --time=48:00:00
#SBATCH --partition=himem

set -euo pipefail

SCRIPT_DIR="/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/VisiumLowDiseaseScripts"
mkdir -p "${SCRIPT_DIR}/logs"
mkdir -p "${SCRIPT_DIR}/output"

module purge
module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0

echo "Hostname  : $(hostname)"
echo "Start time: $(date)"

Rscript "${SCRIPT_DIR}/04_VisiumLowDisease_Deconvolution.R"

echo "Job finished: $(date)"
