#!/bin/bash
#SBATCH --job-name=scRNA_qc
#SBATCH --output=logs/02_qc_%A_%a.out    # logs/ MUST exist before submission (see below)
#SBATCH --error=logs/02_qc_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8     # 8 CPUs: enough for SCTransform + paramSweep
#SBATCH --mem=128G            # doubletFinder(sct=TRUE) runs SCTransform on the
                              # augmented dataset (~30k cells for 24k-cell samples)
                              # and peaks at 20-30 GB; 128G gives safe headroom
#SBATCH --time=8:00:00        # increased: paramSweep on 24k cells takes 2-4 h
#SBATCH --partition=himem
#SBATCH --array=0-22%5       # 23 samples (indices 0–22); %5 = max 5 running at once
                             # Adjust upper bound if samples are added/removed

# IMPORTANT: create logs/ BEFORE submitting this script, not inside it:
#   mkdir -p logs && sbatch submit_02_qc_array.sh
# SLURM opens the log file at job launch, before any commands here run.

module purge
module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# sci_kidney (index 18) has potentially millions of cells — give it more memory
# by detecting its index and re-queuing with higher resources if needed.
# Alternatively, exclude it from this array and submit it separately (recommended):
#   See submit_02_qc_sci.sh

# Pass the array task index to R so each job processes one sample
Rscript 02_qc_and_filter.R --array-index ${SLURM_ARRAY_TASK_ID}
