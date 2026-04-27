#!/bin/bash
#SBATCH --job-name=RenalUro_pipeline
#SBATCH --output=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/RenalUrotheliumScripts/logs/pipeline_%j.out
#SBATCH --error=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/RenalUrotheliumScripts/logs/pipeline_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=1300G
#SBATCH --time=60:00:00
#SBATCH --partition=himem

SCRIPT_DIR="/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/RenalUrotheliumScripts"

mkdir -p "${SCRIPT_DIR}/logs" "${SCRIPT_DIR}/output"

echo "======================================================"
echo "  RenalUrothelium Full Pipeline"
echo "  Host     : $(hostname)"
echo "  Start    : $(date)"
echo "======================================================"
echo ""

# ── Activate Python environment ───────────────────────────────────────────────
source /home/gdbecknelllab/xxw004/Software/miniconda3/etc/profile.d/conda.sh
conda activate cell2loc_env

export SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-16}

# ── Step 1: Load + preprocess all h5ad files ──────────────────────────────────
echo "=== Step 1: Preprocess (load + HVG selection) ==="
echo "Start: $(date)"
python "${SCRIPT_DIR}/01_preprocess.py"
if [ $? -ne 0 ]; then echo "ERROR in Step 1"; exit 1; fi
echo "Done:  $(date)"
echo ""

# ── Step 2: scVI integration ──────────────────────────────────────────────────
echo "=== Step 2: scVI integration ==="
echo "Start: $(date)"
python "${SCRIPT_DIR}/02_scvi_integrate.py"
if [ $? -ne 0 ]; then echo "ERROR in Step 2"; exit 1; fi
echo "Done:  $(date)"
echo ""

# ── Step 3: Extract urothelial cells ─────────────────────────────────────────
echo "=== Step 3: Extract urothelial cells ==="
echo "Start: $(date)"
python "${SCRIPT_DIR}/03_extract_uro.py"
if [ $? -ne 0 ]; then echo "ERROR in Step 3"; exit 1; fi
echo "Done:  $(date)"
echo ""

conda deactivate

# ── Load R environment ────────────────────────────────────────────────────────
module purge
module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0

# Prepend newer libstdc++ so zellkonverter/basilisk finds GLIBCXX_3.4.29
# (GCC/9.3.0 is too old; GCCcore/12.3.0 provides the required version)
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
echo "  Pipeline complete: $(date)"
echo "======================================================"

sacct -j "${SLURM_JOB_ID}" \
  --format=JobID,User,Start,End,Elapsed,TotalCPU,State,NodeList \
  --noheader
