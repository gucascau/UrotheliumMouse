#!/bin/bash
#SBATCH --job-name=bladder_qc
#SBATCH --output=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/BladderUrothelium/logs/02_qc_%A_%a.out
#SBATCH --error=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/BladderUrothelium/logs/02_qc_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=4:00:00
#SBATCH --partition=himem
#SBATCH --array=0-9
# Array indices 0–9 correspond to the 10 bladder samples (sorted alphabetically).
# Run AFTER submit_01_load.sh completes.
#
# To process all samples sequentially instead of as an array:
#   Rscript 02_qc_bladder.R
#
# To regenerate the summary table after array jobs finish:
#   Rscript 02_qc_bladder.R --summary-only

set -euo pipefail

PROJECT_DIR="/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/BladderUrothelium"
SCRIPT_DIR="${PROJECT_DIR}/BladderUrotheliumScripts"

mkdir -p "${PROJECT_DIR}/logs"
mkdir -p "${PROJECT_DIR}/qc_plots"

module purge
module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
: "${SLURM_ARRAY_TASK_ID:?submit_02_qc.sh must be submitted as a SLURM array job}"

echo "Hostname:    $(hostname)"
echo "Start time:  $(date)"
echo "Array index: ${SLURM_ARRAY_TASK_ID}"
echo ""

Rscript "${SCRIPT_DIR}/02_qc_bladder.R" --array-index "${SLURM_ARRAY_TASK_ID}"

echo ""
echo "End time: $(date)"
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  sacct -j "${SLURM_JOB_ID}" \
    --format=JobID,User,Start,End,Elapsed,TotalCPU,State,NodeList \
    --noheader --units=G || true
fi
