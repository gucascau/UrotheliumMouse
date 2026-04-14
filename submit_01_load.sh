#!/bin/bash
#SBATCH --job-name=scRNA_load
#SBATCH --output=logs/01_load_%j.out
#SBATCH --error=logs/01_load_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G          # Loading + sparse matrices for all datasets
#SBATCH --time=12:00:00
#SBATCH --partition=himem

mkdir -p logs

module purge
module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0

# Limit R to use available cores for data.table/Matrix operations
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}

Rscript 01_load_datasets.R
