#!/bin/bash
#SBATCH --job-name=scvi_qc
#SBATCH --output=logs/06_scvi_qc_%j.out
#SBATCH --error=logs/06_scvi_qc_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=384G
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00
#SBATCH --partition=gpu
# Run AFTER 05_export_qc_h5ad.R has written all files in qc_h5ad/

mkdir -p logs

module purge
module load CUDA/12.4.0

source activate cell2loc_env

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

pip install --quiet scikit-misc

python 06_scvi_qc_integration.py
