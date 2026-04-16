#!/bin/bash
#SBATCH --job-name=export_qc_h5ad
#SBATCH --output=logs/05_export_qc_h5ad_%j.out
#SBATCH --error=logs/05_export_qc_h5ad_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G          # EmbryosE9_5ToE13_5 (~1.9 M cells) is the largest sample
#SBATCH --time=06:00:00
#SBATCH --partition=himem

mkdir -p logs

source /home/gdbecknelllab/xxw004/miniconda3/etc/profile.d/conda.sh
conda activate scrna_env

Rscript 05_export_qc_h5ad.R
