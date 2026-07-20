#!/bin/bash
#SBATCH --job-name=urodev_cytotrace
#SBATCH --output=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/UrotheliumDevelopmentScripts/logs/urodev_cytotrace_%j.out
#SBATCH --error=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/UrotheliumDevelopmentScripts/logs/urodev_cytotrace_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --partition=himem

set -euo pipefail

SCRIPT_DIR="/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/UrotheliumDevelopmentScripts"
mkdir -p "${SCRIPT_DIR}/logs"
mkdir -p "${SCRIPT_DIR}/output"

# This script needs slingshot/CytoTRACE, which are only installed in the
# scrna_env conda environment (not the module-loaded R/4.4.0 used by
# 01/02) -- see UrotheliumDevelopmentScripts memory notes.
module purge
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scrna_env

echo "Hostname  : $(hostname)"
echo "Start time: $(date)"

Rscript "${SCRIPT_DIR}/04_UrotheliumDevelopment_CytoTRACE.R"

echo "Job finished: $(date)"
