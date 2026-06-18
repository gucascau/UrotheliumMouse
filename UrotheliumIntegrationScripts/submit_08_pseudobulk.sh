#!/bin/bash
################################################################################
# submit_08_pseudobulk.sh
#
# Pseudobulk DESeq2: compare FinalConditionL2 groups within each
# FinalConditionL1 using sample_id as biological replicates.
#
# Submit with:  sbatch submit_08_pseudobulk.sh
################################################################################

#SBATCH --job-name=AllUro_pseudobulk
#SBATCH --output=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/FinalUrotheliumCells/UrotheliumScripts/logs/AllUro_pseudobulk_%j.out
#SBATCH --error=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/FinalUrotheliumCells/UrotheliumScripts/logs/AllUro_pseudobulk_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=4:00:00
#SBATCH --partition=himem

SCR_DIR="/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/FinalUrotheliumCells/UrotheliumScripts"
mkdir -p "${SCR_DIR}/logs"

module purge
module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0

echo "Hostname  : $(hostname)"
echo "Start time: $(date)"
echo ""

Rscript "${SCR_DIR}/08_pseudobulk_DESeq2.R"

echo ""
echo "End time: $(date)"
sacct -j "${SLURM_JOB_ID}" \
  --format=JobID,User,Start,End,Elapsed,TotalCPU,State,NodeList \
  --noheader --units=G
