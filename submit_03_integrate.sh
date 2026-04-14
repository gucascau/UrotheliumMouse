#!/bin/bash
#SBATCH --job-name=scRNA_integrate
#SBATCH --output=logs/03_integrate_%j.out
#SBATCH --error=logs/03_integrate_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=512G          # Full merge can be large; BPCells on-disk helps
#SBATCH --time=24:00:00
#SBATCH --partition=himem
# No dependency: submit directly once all _qc.rds files are confirmed present.
# (Dependency was: afterok:${QC_JOB_ID} — removed for manual submission)

mkdir -p logs

module purge
module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0

# Set Harmony and Matrix to use all available cores
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}

Rscript 03_integrate_harmony.R
