#!/bin/bash
#SBATCH --job-name=bladder_load
#SBATCH --output=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/BladderUrothelium/logs/01_load_%j.out
#SBATCH --error=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/BladderUrothelium/logs/01_load_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=6:00:00
#SBATCH --partition=himem
# Memory note: GSE164557 files are dense ASCII matrices (200–580 MB each).
# fread loads them fully into RAM before sparsifying; 256 GB covers peak usage.

set -euo pipefail

PROJECT_DIR="/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/BladderUrothelium"
SCRIPT_DIR="${PROJECT_DIR}/BladderUrotheliumScripts"

mkdir -p "${PROJECT_DIR}/logs"
mkdir -p "${PROJECT_DIR}/seurat_objects"

source /home/gdbecknelllab/xxw004/miniconda3/etc/profile.d/conda.sh
conda activate scrna_env

echo "Hostname:   $(hostname)"
echo "Start time: $(date)"
echo ""

Rscript "${SCRIPT_DIR}/01_load_bladder.R"

echo ""
echo "End time: $(date)"
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  sacct -j "${SLURM_JOB_ID}" \
    --format=JobID,User,Start,End,Elapsed,TotalCPU,State,NodeList \
    --noheader --units=G || true
fi
