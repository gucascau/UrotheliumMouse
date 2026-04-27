#!/bin/bash
#SBATCH --job-name=Uro_meta_corr
#SBATCH --output=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/logs/Uro_meta_corr_%j.out
#SBATCH --error=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/logs/Uro_meta_corr_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=2:00:00
#SBATCH --partition=himem

SCRIPT_DIR="/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"

mkdir -p "${SCRIPT_DIR}/logs"

module purge
module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0

echo "Hostname:   $(hostname)"
echo "Start time: $(date)"
echo ""

Rscript "${SCRIPT_DIR}/plot_Urothelium_metadata_correlations.R"

echo ""
echo "End time: $(date)"

sacct -j "${SLURM_JOB_ID}" \
  --format=JobID,User,Start,End,Elapsed,TotalCPU,State,NodeList \
  --noheader --units=G
