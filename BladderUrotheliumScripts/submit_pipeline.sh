#!/bin/bash
#SBATCH --job-name=bladder_pipeline
#SBATCH --output=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/BladderUrothelium/logs/pipeline_%j.out
#SBATCH --error=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/BladderUrothelium/logs/pipeline_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00
#SBATCH --partition=himem

set -euo pipefail

# Submit the bladder urothelium processing pipeline with SLURM dependencies.
#
# Pipeline:
#   01_load -> 02_qc array -> 03_export_h5ad -> 05_scvi
#                        \-> 04_harmony
#
# Usage:
#   sbatch submit_pipeline.sh
#   bash submit_pipeline.sh
#   bash submit_pipeline.sh --dry-run

SCRIPT_DIR="/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/BladderUrothelium/BladderUrotheliumScripts"
PROJECT_DIR="/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/BladderUrothelium"
LOG_DIR="${SCRIPT_DIR}/logs"
PROJECT_LOG_DIR="${PROJECT_DIR}/logs"
DRY_RUN=0

usage() {
  cat <<USAGE
Usage: $(basename "$0") [--dry-run] [--help]

Submits:
  1. submit_01_load.sh
  2. submit_02_qc.sh after successful load
  3. submit_03_export_h5ad.sh after successful QC array
  4. submit_04_harmony.sh after successful QC array
  5. submit_05_scvi.sh after successful h5ad export

Options:
  --dry-run   Print the sbatch commands without submitting jobs.
  --help      Show this help message.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

mkdir -p "${LOG_DIR}" "${PROJECT_LOG_DIR}"

submit_job() {
  local label="$1"
  shift

  if [[ "${DRY_RUN}" -eq 1 ]]; then
    printf '[dry-run] %-12s sbatch' "${label}" >&2
    printf ' %q' "$@" >&2
    printf '\n' >&2
    echo "DRYRUN_${label}"
    return 0
  fi

  local raw_job_id
  raw_job_id="$(sbatch --parsable "$@")"
  local job_id="${raw_job_id%%;*}"
  printf '%s\n' "${job_id}"
}

dependency_arg() {
  local job_id="$1"

  printf '%s\n' "--dependency=afterok:${job_id}"
}

cd "${SCRIPT_DIR}"

echo "Submitting bladder urothelium pipeline from:"
echo "  ${SCRIPT_DIR}"
echo "======================================================"
echo "  BladderUrothelium Pipeline Launcher"
echo "  Host     : $(hostname)"
echo "  Start    : $(date)"
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  if [[ "${DRY_RUN}" -eq 0 ]]; then
    echo "  Launcher : ${SLURM_JOB_ID}"
  fi
fi
echo "======================================================"
echo ""

load_job="$(submit_job load submit_01_load.sh)"
qc_job="$(submit_job qc "$(dependency_arg "${load_job}")" submit_02_qc.sh)"
h5ad_job="$(submit_job h5ad "$(dependency_arg "${qc_job}")" submit_03_export_h5ad.sh)"
harmony_job="$(submit_job harmony "$(dependency_arg "${qc_job}")" submit_04_harmony.sh)"
scvi_job="$(submit_job scvi "$(dependency_arg "${h5ad_job}")" submit_05_scvi.sh)"

echo ""
echo "Submission summary"
printf '%-12s %-18s %s\n' "Step" "Job ID" "Dependency"
printf '%-12s %-18s %s\n' "01_load" "${load_job}" "-"
printf '%-12s %-18s %s\n' "02_qc" "${qc_job}" "afterok:${load_job}"
printf '%-12s %-18s %s\n' "03_h5ad" "${h5ad_job}" "afterok:${qc_job}"
printf '%-12s %-18s %s\n' "04_harmony" "${harmony_job}" "afterok:${qc_job}"
printf '%-12s %-18s %s\n' "05_scvi" "${scvi_job}" "afterok:${h5ad_job}"
echo ""
echo "Check status with:"
echo "  squeue -u ${USER:-$(id -un 2>/dev/null || echo unknown)}"

if [[ -n "${SLURM_JOB_ID:-}" && "${DRY_RUN}" -eq 0 ]]; then
  echo ""
  echo "Launcher job accounting:"
  sacct -j "${SLURM_JOB_ID}" \
    --format=JobID,User,Start,End,Elapsed,TotalCPU,State,NodeList \
    --noheader --units=G || true
fi

echo ""
echo "Pipeline submission complete: $(date)"
