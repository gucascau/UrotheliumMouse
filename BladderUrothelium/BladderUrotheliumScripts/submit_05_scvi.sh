#!/bin/bash
#SBATCH --job-name=bladder_scvi
#SBATCH --output=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/BladderUrothelium/logs/05_scvi_%j.out
#SBATCH --error=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/BladderUrothelium/logs/05_scvi_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --gres=gpu:1
#SBATCH --time=12:00:00
#SBATCH --partition=gpu
# Run AFTER submit_03_export_h5ad.sh completes (reads *.h5ad from qc_h5ad/).

set -euo pipefail

PROJECT_DIR="/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/BladderUrothelium"
SCRIPT_DIR="${PROJECT_DIR}/BladderUrotheliumScripts"

mkdir -p "${PROJECT_DIR}/logs"
mkdir -p "${PROJECT_DIR}/integration_output"

if [[ -n "${LMOD_CMD:-}" ]]; then
  eval "$("${LMOD_CMD}" bash purge)"
  eval "$("${LMOD_CMD}" bash load CUDA/12.4.0)"
else
  module purge
  module load CUDA/12.4.0
fi

source /home/gdbecknelllab/xxw004/Software/miniconda3/etc/profile.d/conda.sh
conda activate cell2loc_env

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export PIPELINE_CACHE_DIR="${TMPDIR:-/tmp}/bladder_scvi_${USER:-user}_${SLURM_JOB_ID:-$$}"
export MPLCONFIGDIR="${PIPELINE_CACHE_DIR}/matplotlib"
export NUMBA_CACHE_DIR="${PIPELINE_CACHE_DIR}/numba"
mkdir -p "${MPLCONFIGDIR}" "${NUMBA_CACHE_DIR}"

python - <<'PY'
import importlib.util
import sys

if importlib.util.find_spec("skmisc") is None:
    sys.exit("ERROR: scikit-misc is required in cell2loc_env for Scanpy HVG selection.")
PY

echo "Hostname:   $(hostname)"
echo "Start time: $(date)"
echo ""

python "${SCRIPT_DIR}/05_scvi_bladder.py"

echo ""
echo "End time: $(date)"
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  sacct -j "${SLURM_JOB_ID}" \
    --format=JobID,User,Start,End,Elapsed,TotalCPU,State,NodeList \
    --noheader --units=G || true
fi
