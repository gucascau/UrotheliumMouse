#!/bin/bash
#SBATCH --job-name=urodev_bigwig
#SBATCH --output=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/UrotheliumDevelopmentScripts/logs/urodev_bigwig_%j.out
#SBATCH --error=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/UrotheliumDevelopmentScripts/logs/urodev_bigwig_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=160G
#SBATCH --time=06:00:00
#SBATCH --partition=himem

set -euo pipefail

SCRIPT_DIR="/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/UrotheliumDevelopmentScripts"
mkdir -p "${SCRIPT_DIR}/logs"
mkdir -p "${SCRIPT_DIR}/output"

module purge
module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0

echo "Hostname  : $(hostname)"
echo "Start time: $(date)"

Rscript "${SCRIPT_DIR}/13_UrotheliumDevelopment_GenePanel_BigWig.R"

echo "Job finished: $(date)"
