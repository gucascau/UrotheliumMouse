#!/bin/bash
#SBATCH --job-name=export_merged_h5ad
#SBATCH --output=logs/03c_export_merged_%j.out
#SBATCH --error=logs/03c_export_merged_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=512G          # merged_normalized.rds is ~52 GB on disk; JoinLayers + transpose needs ~400-512 GB RAM
#SBATCH --time=06:00:00
#SBATCH --partition=himem
# Run AFTER 03_integrate_harmony.R has written:
#   integration_output/merged_normalized.rds
#   integration_output/hvg_list.rds

mkdir -p logs

source /home/gdbecknelllab/xxw004/miniconda3/etc/profile.d/conda.sh
conda activate scrna_env
# module purge
# module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0

# # Point reticulate to the cell2loc_env Python (has anndata installed)
# export RETICULATE_PYTHON="/home/gdbecknelllab/xxw004/Workspace/.conda/envs/cell2loc_env/bin/python"

# # Use the conda env's libstdc++.so.6 instead of GCC/9.3.0's (which lacks GLIBCXX_3.4.29).
# # pandas >=1.4 requires GLIBCXX_3.4.29, but GCCcore/9.3.0 only ships up to GLIBCXX_3.4.28.
# export LD_LIBRARY_PATH="/home/gdbecknelllab/xxw004/Workspace/.conda/envs/cell2loc_env/lib:${LD_LIBRARY_PATH}"



Rscript 03c_export_merged_normalized_h5ad.R
