#!/bin/bash
################################################################################
# submit_04_scvi.sh
#
# Two-step job:
#   Step 1 (R)      : export gated Seurat RDS to sparse MTX + metadata
#   Step 2 (Python) : train scVI, compute UMAP + Leiden, save h5ad
#
# Submit with:  sbatch submit_04_scvi.sh
#
# To use a GPU node (much faster training — ~15 min vs ~3 h), change
# --partition and add --gres=gpu:1 then uncomment the GPU env lines below.
################################################################################

#SBATCH --job-name=AllUro_scvi
#SBATCH --output=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/FinalUrotheliumCells/UrotheliumScripts/logs/AllUro_scvi_%j.out
#SBATCH --error=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/FinalUrotheliumCells/UrotheliumScripts/logs/AllUro_scvi_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=300G
#SBATCH --time=12:00:00
#SBATCH --partition=himem
# To use GPU instead:
#   --partition=gpu  --gres=gpu:1  --mem=128G  --time=4:00:00

SCR_DIR="/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/FinalUrotheliumCells/UrotheliumScripts"
mkdir -p "${SCR_DIR}/logs"

echo "Hostname  : $(hostname)"
echo "Start time: $(date)"
echo ""

# # ── Step 1: R export ───────────────────────────────────────────────────────────
# echo "===== Step 1: Export RDS to scvi_input ====="
# module purge
# module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0

# Rscript "${SCR_DIR}/04a_export_for_scvi.R"
# R_EXIT=$?

# if [ ${R_EXIT} -ne 0 ]; then
#   echo "ERROR: R export failed (exit ${R_EXIT}) — aborting." >&2
#   exit ${R_EXIT}
# fi

echo ""
echo "===== Step 2: scVI training ====="

# ── Step 2: Python scVI ────────────────────────────────────────────────────────
# Use the full path to Python in the conda env — conda activate does not
# reliably modify the environment in non-interactive SLURM bash sessions.
PYTHON="/home/gdbecknelllab/xxw004/Workspace/.conda/envs/cell2loc_env/bin/python3"

# Set number of threads for PyTorch (CPU training)
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# Uncomment for GPU runs to restrict to the allocated GPU:
# export CUDA_VISIBLE_DEVICES=$(nvidia-smi --query-gpu=index --format=csv,noheader | head -1)

"${PYTHON}" "${SCR_DIR}/04b_scvi_AllUrothelium.py"
PY_EXIT=$?

echo ""
echo "End time: $(date)"
sacct -j "${SLURM_JOB_ID}" \
  --format=JobID,User,Start,End,Elapsed,TotalCPU,State,NodeList \
  --noheader --units=G

exit ${PY_EXIT}
