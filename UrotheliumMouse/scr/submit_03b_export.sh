#!/bin/bash
#SBATCH --job-name=scVI_export
#SBATCH --output=logs/03b_export_%j.out
#SBATCH --error=logs/03b_export_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=512G          # merged_pre_integration.rds is 54 GB compressed; needs ~400 GB in R
#SBATCH --time=06:00:00
#SBATCH --partition=himem
# Run AFTER 03_integrate_harmony.R has written merged_pre_integration.rds

mkdir -p logs

module purge
module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0

# Point reticulate to the cell2loc_env Python (has anndata installed)
export RETICULATE_PYTHON="/home/gdbecknelllab/xxw004/Workspace/.conda/envs/cell2loc_env/bin/python"

Rscript 03b_export_for_scvi.R
