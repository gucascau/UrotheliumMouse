#!/bin/bash
#SBATCH --job-name=scVI_export
#SBATCH --output=logs/03b_export_%j.out
#SBATCH --error=logs/03b_export_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G          # reads hvg_list.rds (tiny) + per-sample QC RDS files; no need to load full merged object
#SBATCH --time=06:00:00
#SBATCH --partition=himem
# Run AFTER 03_integrate_harmony.R has written hvg_list.rds

mkdir -p logs

module purge
module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0

# Point reticulate to the cell2loc_env Python (has anndata installed)
export RETICULATE_PYTHON="/home/gdbecknelllab/xxw004/Workspace/.conda/envs/cell2loc_env/bin/python"

# Use the conda env's libstdc++.so.6 instead of GCC/9.3.0's (which lacks GLIBCXX_3.4.29).
# pandas >=1.4 requires GLIBCXX_3.4.29, but GCCcore/9.3.0 only ships up to GLIBCXX_3.4.28.
export LD_LIBRARY_PATH="/home/gdbecknelllab/xxw004/Workspace/.conda/envs/cell2loc_env/lib:${LD_LIBRARY_PATH}"

Rscript 03b_export_for_scvi.R
