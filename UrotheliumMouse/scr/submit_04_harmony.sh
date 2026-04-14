#!/bin/bash
#SBATCH --job-name=harmony_integrate
#SBATCH --output=logs/04_harmony_%j.out
#SBATCH --error=logs/04_harmony_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=512G          # merged_normalized.rds is large; ScaleData + PCA need ~400 GB
#SBATCH --time=12:00:00
#SBATCH --partition=himem
# Run AFTER 03_integrate_harmony.R has written:
#   integration_output/merged_normalized.rds
#   integration_output/hvg_list.rds

mkdir -p logs

module purge
module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0

Rscript 04_integrate_harmony.R
