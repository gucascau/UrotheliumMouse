#!/bin/bash
#SBATCH --job-name=fix_lapack
#SBATCH --output=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/FinalUrotheliumCells/UrotheliumScripts/logs/fix_lapack_%j.out
#SBATCH --error=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/FinalUrotheliumCells/UrotheliumScripts/logs/fix_lapack_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00

SCR_DIR="/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/FinalUrotheliumCells/UrotheliumScripts"
mkdir -p "${SCR_DIR}/logs"

module purge
module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0

echo "Hostname  : $(hostname)"
echo "Start time: $(date)"
echo ""

Rscript "${SCR_DIR}/fix_lapack_packages.R"

echo ""
echo "End time: $(date)"
echo ""
echo "If the fix succeeded, resubmit with:"
echo "  sbatch ${SCR_DIR}/submit_01_integrate.sh"
