#!/bin/bash
#SBATCH --job-name=probe_rctd_purity
#SBATCH --output=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/XeniumDevelopmentScripts/logs/probe_rctd_purity_%j.out
#SBATCH --error=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/XeniumDevelopmentScripts/logs/probe_rctd_purity_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=150G
#SBATCH --time=00:30:00
#SBATCH --partition=himem

module purge
module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0
Rscript "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/XeniumDevelopmentScripts/probe_rctd_purity.R"
