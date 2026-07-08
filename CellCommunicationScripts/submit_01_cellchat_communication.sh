#!/bin/bash
#SBATCH --job-name=Uro_CellChat
#SBATCH --output=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/CellCommunicationScripts/logs/CellChat_Urothelium_%j.out
#SBATCH --error=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/CellCommunicationScripts/logs/CellChat_Urothelium_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=1500G
#SBATCH --time=72:00:00
#SBATCH --partition=himem
# Full-scale CellChat run (no downsampling) on the 968k-cell
# RenalUrothelium_allcells_scvi_annotations_meta_uro_updated.rds object,
# focused on Urothelium <-> other kidney cell type communication.

SCR_DIR="/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/CellCommunicationScripts"
mkdir -p "${SCR_DIR}/logs"
mkdir -p "${SCR_DIR}/output"

source /home/gdbecknelllab/xxw004/miniconda3/etc/profile.d/conda.sh
conda activate scrna_env

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo "Hostname:   $(hostname)"
echo "Start time: $(date)"
echo ""

Rscript "${SCR_DIR}/01_Urothelium_CellChat_Communication.R"

echo ""
echo "End time: $(date)"
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  sacct -j "${SLURM_JOB_ID}" \
    --format=JobID,User,Start,End,Elapsed,TotalCPU,State,NodeList \
    --noheader --units=G || true
fi
