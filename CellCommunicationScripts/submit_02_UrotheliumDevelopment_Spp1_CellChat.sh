#!/bin/bash
#SBATCH --job-name=UroDev_Spp1_CellChat
#SBATCH --output=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/CellCommunicationScripts/logs/UroDev_Spp1_CellChat_%j.out
#SBATCH --error=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/CellCommunicationScripts/logs/UroDev_Spp1_CellChat_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --partition=himem
# 6 sequential per-stage CellChat runs (Chen2025 developmental atlas,
# ~34k cells/stage on average, full celltype_final grouping) -- much smaller
# per-run than submit_01_cellchat_communication.sh's single 968k-cell run,
# so a fraction of its resource request, but 6 runs in series.

SCR_DIR="/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/CellCommunicationScripts"
mkdir -p "${SCR_DIR}/logs"
mkdir -p "${SCR_DIR}/output"

source /home/gdbecknelllab/xxw004/miniconda3/etc/profile.d/conda.sh
conda activate scrna_env

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo "Hostname:   $(hostname)"
echo "Start time: $(date)"
echo ""

Rscript "${SCR_DIR}/02_UrotheliumDevelopment_Spp1_CellChat.R"

echo ""
echo "End time: $(date)"
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  sacct -j "${SLURM_JOB_ID}" \
    --format=JobID,User,Start,End,Elapsed,TotalCPU,State,NodeList \
    --noheader --units=G || true
fi
