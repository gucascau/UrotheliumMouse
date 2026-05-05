#!/bin/bash
#SBATCH --job-name=bladder_harmony
#SBATCH --output=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/BladderUrothelium/logs/04_harmony_%j.out
#SBATCH --error=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/BladderUrothelium/logs/04_harmony_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=1170G
#SBATCH --time=24:00:00
#SBATCH --partition=himem
# Run AFTER submit_02_qc.sh completes (reads *_qc.rds files).
# Memory: ScaleData on 3000 HVGs across ~100k cells needs ~64 GB;
#         256 GB provides comfortable headroom.

set -euo pipefail

PROJECT_DIR="/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/BladderUrothelium"
SCRIPT_DIR="${PROJECT_DIR}/BladderUrotheliumScripts"

mkdir -p "${PROJECT_DIR}/logs"
mkdir -p "${PROJECT_DIR}/integration_output"

source /home/gdbecknelllab/xxw004/miniconda3/etc/profile.d/conda.sh
conda activate scrna_env

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo "Hostname:   $(hostname)"
echo "Start time: $(date)"
echo ""

Rscript "${SCRIPT_DIR}/04_harmony_bladder.R"

echo ""
echo "End time: $(date)"
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  sacct -j "${SLURM_JOB_ID}" \
    --format=JobID,User,Start,End,Elapsed,TotalCPU,State,NodeList \
    --noheader --units=G || true
fi
