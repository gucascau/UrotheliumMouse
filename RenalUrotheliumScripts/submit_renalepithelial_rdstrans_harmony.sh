#!/bin/bash
#SBATCH --job-name=RenalEpithelial_Transfer
#SBATCH --output=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/RenalUrotheliumScripts/logs/R_RenalEpithelial_%j.out
#SBATCH --error=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/RenalUrotheliumScripts/logs/R_RenalEpithelial_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=1340G
#SBATCH --time=48:00:00
#SBATCH --partition=himem


SCRIPT_DIR="/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/RenalUrotheliumScripts"

mkdir -p "${SCRIPT_DIR}/logs"


# ── Load R environment ────────────────────────────────────────────────────────
module purge
module load GCC/9.3.0 OpenMPI/4.0.3 OpenBLAS/0.3.9 R/4.4.0

# Prepend newer libstdc++ so zellkonverter/basilisk finds GLIBCXX_3.4.29
export LD_LIBRARY_PATH=/gpfs0/export/apps/easybuild/software/GCCcore/12.3.0/lib64:${LD_LIBRARY_PATH}

echo "======================================================"
echo "  RenalUrothelium Transfer + R Steps "
echo "  Host     : $(hostname)"
echo "  Start    : $(date)"
echo "======================================================"
echo ""

Rscript "${SCRIPT_DIR}/04b_convert_renal_epithelial_to_rds.R"

echo "Renal Epithelial Transfer + R Steps Completed!"
echo "End      : $(date)"
echo "======================================================"

# ── Step 5: Harmony reclustering ─────────────────────────────────────────────

echo "======================================================"
echo "  Step 5: Harmony reclustering "
echo "  Host     : $(hostname)"
echo "  Start    : $(date)"
echo "======================================================"
echo ""

Rscript "${SCRIPT_DIR}/05b_harmony_recluster_renal_epithelial.R"

echo "Renal Epithelial Transfer + R Steps Completed!"
echo "End      : $(date)"
echo "======================================================"