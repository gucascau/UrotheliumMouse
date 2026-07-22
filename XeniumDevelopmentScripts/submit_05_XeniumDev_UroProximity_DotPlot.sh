#!/bin/bash
#SBATCH --job-name=xeniumdev_uroprox
#SBATCH --output=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/XeniumDevelopmentScripts/logs/xeniumdev_uroprox_%j.out
#SBATCH --error=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/XeniumDevelopmentScripts/logs/xeniumdev_uroprox_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=150G
#SBATCH --time=02:00:00
#SBATCH --partition=himem

set -euo pipefail

SCRIPT_DIR="/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/XeniumDevelopmentScripts"
mkdir -p "${SCRIPT_DIR}/logs"
mkdir -p "${SCRIPT_DIR}/output"

module purge
module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0

echo "Hostname  : $(hostname)"
echo "Start time: $(date)"

Rscript "${SCRIPT_DIR}/05_XeniumDev_UroProximity_DotPlot.R"

echo "Job finished: $(date)"
