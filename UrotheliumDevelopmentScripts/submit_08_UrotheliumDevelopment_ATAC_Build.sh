#!/bin/bash
#SBATCH --job-name=urodev_atac_build
#SBATCH --output=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/UrotheliumDevelopmentScripts/logs/urodev_atac_build_%j.out
#SBATCH --error=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/UrotheliumDevelopmentScripts/logs/urodev_atac_build_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --partition=himem

set -euo pipefail

SCRIPT_DIR="/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/UrotheliumDevelopmentScripts"
mkdir -p "${SCRIPT_DIR}/logs"
mkdir -p "${SCRIPT_DIR}/output"

# Hybrid module+conda: Signac/EnsDb.Mmusculus.v79/BSgenome.Mmusculus.UCSC.mm10
# live under module R/4.4.0, but the MACS3 binary Signac::CallPeaks() shells
# out to only exists in the scrna_env conda env. scrna_env ships its OWN
# R/Rscript -- prepending (or even appending) its bin/ to PATH risks the
# wrong R getting picked up by anything that does a bare `Rscript`/`R` PATH
# lookup, so conda is never activated and its bin/ never touches PATH here.
# Instead MACS3_BIN is set to the conda env's binary by absolute path, and
# 08_UrotheliumDevelopment_ATAC_Build.R reads that directly (Sys.getenv)
# rather than relying on `which macs3` -- caught interactively (a smoke
# test of this exact PATH-prepend approach silently broke Rscript itself)
# before it ever reached a submitted job.
module purge
module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0
export MACS3_BIN="/home/gdbecknelllab/xxw004/Workspace/.conda/envs/scrna_env/bin/macs3"

echo "Hostname  : $(hostname)"
echo "Start time: $(date)"
echo "macs3     : ${MACS3_BIN}"

Rscript "${SCRIPT_DIR}/08_UrotheliumDevelopment_ATAC_Build.R"

echo "Job finished: $(date)"
