#!/bin/bash
#SBATCH --job-name=scRNA_pipeline
#SBATCH --output=logs/pipeline_%j.out
#SBATCH --error=logs/pipeline_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G              # this script only submits jobs, no computation
#SBATCH --time=0:10:00        # just needs a few seconds to submit the chain
#SBATCH --partition=himem

# submit_pipeline.sh
# Submit the full integration pipeline as a SLURM dependency chain.
# Usage: sbatch submit_pipeline.sh
#
# Pipeline:
#   Step 1 – Load all datasets           (256 GB, 8 CPUs, 12 h)
#   Step 2a – QC array (all except sci)  (64 GB × 5 parallel, 8 CPUs, 4 h)
#   Step 2b – QC sci_kidney separately   (256 GB, 8 CPUs, 12 h)
#   Step 2c – Write summary table        (sequential, after 2a+2b)
#   Step 3 – Harmony integration         (512 GB, 16 CPUs, 24 h)

set -euo pipefail

# logs/ must exist before SLURM can write output files
mkdir -p logs

# ── Step 1: Load datasets – SKIPPED (already completed successfully) ──────────
# JOB1=$(sbatch --parsable submit_01_load.sh)
# echo "Submitted Step 1 (load): job ${JOB1}"

# ── Step 2a: QC array – all samples except EmbryosE9_5ToE13_5 ─────────────────
# Compute the 0-based index of EmbryosE9_5ToE13_5 at runtime so this stays
# correct even if samples are added/removed (do NOT hardcode).
SCI_IDX=$(Rscript --vanilla -e '
  OBJ_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/seurat_objects"
  f <- sort(list.files(OBJ_DIR, pattern = "_seurat\\.rds$"))
  cat(which(sub("_seurat\\.rds$", "", f) == "EmbryosE9_5ToE13_5") - 1L)
')
N_MAX=$(Rscript --vanilla -e '
  OBJ_DIR <- "/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/seurat_objects"
  cat(length(list.files(OBJ_DIR, pattern = "_seurat\\.rds$")) - 1L)
')
# Build array range excluding SCI_IDX
if   [ "${SCI_IDX}" -eq 0 ]; then
  ARRAY_STR="1-${N_MAX}"
elif [ "${SCI_IDX}" -eq "${N_MAX}" ]; then
  ARRAY_STR="0-$((N_MAX - 1))"
else
  ARRAY_STR="0-$((SCI_IDX - 1)),$((SCI_IDX + 1))-${N_MAX}"
fi
echo "EmbryosE9_5ToE13_5 is at index ${SCI_IDX}; running array: ${ARRAY_STR}"

JOB2A=$(sbatch --parsable \
               --array=${ARRAY_STR}%5 \
               submit_02_qc_array.sh)
echo "Submitted Step 2a (QC array ${ARRAY_STR}): job ${JOB2A}"

# ── Step 2b: QC sci_kidney (256 GB → himem) ───────────────────────────────────
JOB2B=$(sbatch --parsable \
               submit_02_qc_sci.sh)
echo "Submitted Step 2b (QC sci_kidney):    job ${JOB2B}"

# ── Step 2c: Write QC summary table after 2a + 2b finish (general) ────────────
JOB2C=$(sbatch --dependency=afterok:${JOB2A}:${JOB2B} \
               --parsable \
               --partition=himem \
               --job-name=scRNA_qc_summary \
               --output=logs/02_qc_summary_%j.out \
               --error=logs/02_qc_summary_%j.err \
               --ntasks=1 --cpus-per-task=2 --mem=32G --time=1:00:00 \
               --wrap="module purge && module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0 && Rscript 02_qc_and_filter.R --summary-only")
echo "Submitted Step 2c (QC summary):       job ${JOB2C} (depends on ${JOB2A} + ${JOB2B})"

# ── Step 3: Harmony integration (512 GB → himem) ─────────────────────────────
JOB3=$(sbatch --dependency=afterok:${JOB2C} \
              --parsable \
              submit_03_integrate.sh)
echo "Submitted Step 3 (integrate):         job ${JOB3} (depends on ${JOB2C})"

echo ""
echo "Full pipeline submitted. Monitor with: squeue -u \$USER"
echo "Logs in: logs/"
echo ""
echo "Job IDs:"
echo "  Step 2a (QC array):    ${JOB2A}"
echo "  Step 2b (QC sci):      ${JOB2B}"
echo "  Step 2c (QC summary):  ${JOB2C}"
echo "  Step 3 (integrate):    ${JOB3}"
