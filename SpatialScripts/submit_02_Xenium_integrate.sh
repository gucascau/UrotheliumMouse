#!/bin/bash
#SBATCH --job-name=xenium_integrate
#SBATCH --output=logs/xenium_integrate_%j.out
#SBATCH --error=logs/xenium_integrate_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=700G
#SBATCH --time=24:00:00
#SBATCH --partition=himem

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "${SCRIPT_DIR}/logs"
mkdir -p "${SCRIPT_DIR}/output"

echo "Job started: $(date)"
echo "Node: ${SLURM_NODELIST}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"

module purge
module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0

Rscript "${SCRIPT_DIR}/02_Xenium_integrate_harmony.R"

echo "Job finished: $(date)"
