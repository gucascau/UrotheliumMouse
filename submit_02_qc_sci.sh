#!/bin/bash
#SBATCH --job-name=scRNA_qc_sci
#SBATCH --output=logs/02_qc_sci_%j.out
#SBATCH --error=logs/02_qc_sci_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G            # EmbryosE9_5ToE13_5 has ~1.9M cells; loading alone needs himem
#SBATCH --time=12:00:00
#SBATCH --partition=himem

mkdir -p logs

module purge
module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# Resolve sci_kidney's 0-based array index at runtime so it stays correct
# if samples are ever added or removed (avoids hardcoding index 18).
Rscript 02_qc_and_filter.R --array-index $(
  Rscript --vanilla -e '
    DATA_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells"
    OBJ_DIR  <- file.path(DATA_DIR, "seurat_objects")
    f <- sort(list.files(OBJ_DIR, pattern = "_seurat\\.rds$"))
    cat(which(sub("_seurat\\.rds$","",f) == "EmbryosE9_5ToE13_5") - 1L)
  ')
