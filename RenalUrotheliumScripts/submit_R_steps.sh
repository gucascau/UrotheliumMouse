#!/bin/bash
#SBATCH --job-name=RenalUro_extract_R
#SBATCH --output=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/RenalUrotheliumScripts/logs/R_steps_%j.out
#SBATCH --error=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/RenalUrotheliumScripts/logs/R_steps_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=1300G
#SBATCH --time=24:00:00
#SBATCH --partition=himem

SCRIPT_DIR="/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/RenalUrotheliumScripts"

mkdir -p "${SCRIPT_DIR}/logs"

echo "======================================================"
echo "  RenalUrothelium Extract + R Steps (3 + 4 + 5)"
echo "  Host     : $(hostname)"
echo "  Start    : $(date)"
echo "======================================================"
echo ""

# ── Step 3: Extract urothelial cells ─────────────────────────────────────────
echo "=== Step 3: Extract urothelial cells ==="
echo "Start: $(date)"

source /home/gdbecknelllab/xxw004/Software/miniconda3/etc/profile.d/conda.sh
conda activate cell2loc_env
export SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-32}

python "${SCRIPT_DIR}/03_extract_uro.py"
if [ $? -ne 0 ]; then echo "ERROR in Step 3"; exit 1; fi

conda deactivate
echo "Done:  $(date)"
echo ""

# ── Load R environment ────────────────────────────────────────────────────────
module purge
module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0

# Prepend newer libstdc++ so zellkonverter/basilisk finds GLIBCXX_3.4.29
export LD_LIBRARY_PATH=/gpfs0/export/apps/easybuild/software/GCCcore/12.3.0/lib64:${LD_LIBRARY_PATH}

# ── Step 4: Convert h5ad → Seurat RDS ────────────────────────────────────────
echo "=== Step 4: Convert h5ad → Seurat RDS ==="
echo "Start: $(date)"
Rscript "${SCRIPT_DIR}/04_convert_to_rds.R"
if [ $? -ne 0 ]; then echo "ERROR in Step 4"; exit 1; fi
echo "Done:  $(date)"
echo ""

# ── Step 5: Harmony reclustering ─────────────────────────────────────────────
echo "=== Step 5: Harmony reclustering ==="
echo "Start: $(date)"
Rscript "${SCRIPT_DIR}/05_harmony_recluster.R"
if [ $? -ne 0 ]; then echo "ERROR in Step 5"; exit 1; fi
echo "Done:  $(date)"
echo ""

echo "======================================================"
echo "  Extract + R steps complete: $(date)"
echo "======================================================"

sacct -j "${SLURM_JOB_ID}" \
  --format=JobID,User,Start,End,Elapsed,TotalCPU,State,NodeList \
  --noheader
