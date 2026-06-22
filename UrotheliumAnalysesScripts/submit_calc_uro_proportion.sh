#!/bin/bash
#
# submit_calc_uro_proportion.sh
#
# Two-step pipeline to compute urothelium proportions across all datasets:
#
#   Step 1 (Python) — calc_uro_proportion_h5ad.py
#     Reads obs metadata from kidney and bladder h5ad files without loading the
#     expression matrix.  Outputs:
#       output/uro_proportion/uro_proportion_kidney.csv
#       output/uro_proportion/uro_proportion_bladder.csv
#
#   Step 2 (R) — calc_uro_proportion_organoids.R
#     Loads small organoid _qc.rds files, applies the urothelium marker filter,
#     infers KudoUUOUrothelium counts from the 1.7 GB integrated object, then
#     combines all results into:
#       output/uro_proportion/uro_proportion_organoids.csv
#       output/uro_proportion/uro_proportion_combined.csv
#       output/uro_proportion/uro_proportion_plots.pdf
#
# Usage:
#   bash submit_calc_uro_proportion.sh        # submit both steps
#   bash submit_calc_uro_proportion.sh py     # submit Python step only
#   bash submit_calc_uro_proportion.sh r      # submit R step only (no dependency)
#
##############################################################################

SCR_DIR="/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/UrotheliumIntegrationScripts"
LOG_DIR="${SCR_DIR}/logs"
mkdir -p "${LOG_DIR}"

MODE="${1:-both}"

# ── Step 1: Python (h5ad obs reading, very low memory) ───────────────────────
if [[ "${MODE}" == "both" || "${MODE}" == "py" ]]; then
  JOB1=$(sbatch --parsable \
    --job-name=uro_prop_h5ad \
    --output="${LOG_DIR}/uro_prop_h5ad_%j.out" \
    --error="${LOG_DIR}/uro_prop_h5ad_%j.err" \
    --ntasks=1 \
    --cpus-per-task=4 \
    --mem=32G \
    --time=1:00:00 \
    --partition=normal \
    --wrap="
      module purge
      module load GCC/9.3.0 OpenMPI/4.0.3
      # activate conda env that has scanpy/h5py/pandas
      source \$(conda info --base)/etc/profile.d/conda.sh
      conda activate scvi-env 2>/dev/null || conda activate scanpy 2>/dev/null || true

      echo 'Hostname  : \$(hostname)'
      echo 'Start time: \$(date)'
      echo ''
      python3 '${SCR_DIR}/calc_uro_proportion_h5ad.py'
      echo ''
      echo 'End time: \$(date)'
      sacct -j \"\${SLURM_JOB_ID}\" \
        --format=JobID,User,Start,End,Elapsed,TotalCPU,State,NodeList \
        --noheader --units=G
    "
  )
  echo "Submitted Python step: job ${JOB1}"
fi

# ── Step 2: R (organoid RDS + combine) ───────────────────────────────────────
# Memory: ~16 GB for small organoids + 1.7 GB AllUrothelium_harmony_integrated.rds
if [[ "${MODE}" == "both" || "${MODE}" == "r" ]]; then
  DEP_ARG=""
  if [[ "${MODE}" == "both" && -n "${JOB1}" ]]; then
    DEP_ARG="--dependency=afterok:${JOB1}"
  fi

  JOB2=$(sbatch --parsable \
    --job-name=uro_prop_organoids \
    --output="${LOG_DIR}/uro_prop_organoids_%j.out" \
    --error="${LOG_DIR}/uro_prop_organoids_%j.err" \
    --ntasks=1 \
    --cpus-per-task=8 \
    --mem=64G \
    --time=2:00:00 \
    --partition=normal \
    ${DEP_ARG} \
    --wrap="
      module purge
      module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0

      echo 'Hostname  : \$(hostname)'
      echo 'Start time: \$(date)'
      echo ''
      Rscript '${SCR_DIR}/calc_uro_proportion_organoids.R'
      echo ''
      echo 'End time: \$(date)'
      sacct -j \"\${SLURM_JOB_ID}\" \
        --format=JobID,User,Start,End,Elapsed,TotalCPU,State,NodeList \
        --noheader --units=G
    "
  )
  echo "Submitted R step    : job ${JOB2}"
fi

echo ""
echo "Outputs will be written to:"
echo "  ${SCR_DIR}/output/uro_proportion/"
echo ""
echo "Note: KudoUUOUrothelium_qc.rds (46 GB) is NOT loaded in the R step."
echo "  Cell counts for KUDO samples are inferred from AllUrothelium_harmony_integrated.rds."
echo "  To apply the marker filter to KudoUUOUrothelium directly, request:"
echo "    --mem=200G --partition=himem"
echo "  and set LOAD_KUDO_FULL=TRUE in calc_uro_proportion_organoids.R."
