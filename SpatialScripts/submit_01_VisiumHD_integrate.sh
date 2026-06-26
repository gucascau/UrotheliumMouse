#!/bin/bash
#SBATCH --job-name=visiumhd_integrate
#SBATCH --output=/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/SpatialScripts/logs/visiumhd_integrate_%j.out
#SBATCH --error=/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/SpatialScripts/logs/visiumhd_integrate_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=700G
#SBATCH --time=18:00:00
#SBATCH --partition=himem

set -euo pipefail

SCRIPT_DIR="/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/SpatialScripts"
mkdir -p "${SCRIPT_DIR}/logs"
mkdir -p "${SCRIPT_DIR}/output"

echo "Job started: $(date)"
echo "Node: ${SLURM_NODELIST}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"

module purge
module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0  # adjust to available R module

Rscript "${SCRIPT_DIR}/01_VisiumHD_integrate_harmony_sketch.R"

echo "Job finished: $(date)"
