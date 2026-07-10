#!/bin/bash
#SBATCH --job-name=visiumhd_fix
#SBATCH --output=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/SpatialScripts/logs/visiumhd_fix_%j.out
#SBATCH --error=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/SpatialScripts/logs/visiumhd_fix_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=2:00:00
#SBATCH --partition=himem

# Recovery script — loads RCTD checkpoint, skips the 1.5-day RCTD compute.
# Fixes: AddMetaData crash when RCTD drops low-UMI bins in "full" mode.

SCRIPT_DIR="/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/SpatialScripts"
mkdir -p "${SCRIPT_DIR}/logs"

module purge
module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0

echo "Hostname  : $(hostname)"
echo "Start time: $(date)"
echo ""

Rscript "${SCRIPT_DIR}/02b_fix_AddMetadata.R"

echo ""
echo "End time: $(date)"
sacct -j "${SLURM_JOB_ID}" \
  --format=JobID,User,Start,End,Elapsed,TotalCPU,State,NodeList \
  --noheader --units=G
