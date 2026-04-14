#!/bin/bash
#SBATCH --job-name=scVI_integrate
#SBATCH --output=logs/04_scvi_%j.out
#SBATCH --error=logs/04_scvi_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --gres=gpu:1            # scVI trains ~10x faster on GPU
#SBATCH --time=12:00:00
#SBATCH --partition=gpu         # adjust to your cluster's GPU partition name
# Run AFTER 03b_export_for_scvi.R has written scvi_input.h5ad

mkdir -p logs

module purge
module load CUDA/12.4.0   # matches torch 2.5.1+cu124 in cell2loc_env

# cell2loc_env has: scvi-tools 1.2.1, scanpy 1.10.4, torch 2.5.1+cu124
source activate cell2loc_env

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

python 04_scvi_integration.py
