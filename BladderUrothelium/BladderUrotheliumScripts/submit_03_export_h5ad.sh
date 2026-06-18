#!/bin/bash
#SBATCH --job-name=bladder_h5ad
#SBATCH --output=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/BladderUrothelium/logs/03_export_h5ad_%j.out
#SBATCH --error=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/BladderUrothelium/logs/03_export_h5ad_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --partition=himem
# Run AFTER submit_02_qc.sh (all array tasks) completes.

set -euo pipefail

PROJECT_DIR="/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/BladderUrothelium"
SCRIPT_DIR="${PROJECT_DIR}/BladderUrotheliumScripts"

mkdir -p "${PROJECT_DIR}/logs"
mkdir -p "${PROJECT_DIR}/qc_h5ad"

source /home/gdbecknelllab/xxw004/miniconda3/etc/profile.d/conda.sh
conda activate scrna_env

echo "Hostname:   $(hostname)"
echo "Start time: $(date)"
echo ""

Rscript "${SCRIPT_DIR}/03_export_h5ad_bladder.R"

echo ""
echo "End time: $(date)"
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  sacct -j "${SLURM_JOB_ID}" \
    --format=JobID,User,Start,End,Elapsed,TotalCPU,State,NodeList \
    --noheader --units=G || true
fi
