#!/bin/bash
#SBATCH --job-name=scRNA_qc
#SBATCH --output=logs/02_qc_%j.out
#SBATCH --error=logs/02_qc_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G           # QC reads one object at a time, lower memory
#SBATCH --time=4:00:00
#SBATCH --partition=your_partition    # <-- change to your cluster partition
#SBATCH --dependency=afterok:${LOAD_JOB_ID}  # set LOAD_JOB_ID from job 01

mkdir -p logs

module purge
module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}

Rscript 02_qc_and_filter.R
