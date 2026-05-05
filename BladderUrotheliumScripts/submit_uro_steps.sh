#!/bin/bash
#SBATCH --job-name=BladderUro_steps
#SBATCH --output=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/BladderUrothelium/BladderUrotheliumScripts/logs/uro_steps_%j.out
#SBATCH --error=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/BladderUrothelium/BladderUrotheliumScripts/logs/uro_steps_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=400G
#SBATCH --time=48:00:00
#SBATCH --partition=himem

SCRIPT_DIR="/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/BladderUrothelium/BladderUrotheliumScripts"

mkdir -p "${SCRIPT_DIR}/logs" "${SCRIPT_DIR}/output"

echo "======================================================"
echo "  BladderUrothelium Full-Gene scVI + Extract + R Steps"
echo "  Host     : $(hostname)"
echo "  Start    : $(date)"
echo "======================================================"
echo ""

source /home/gdbecknelllab/xxw004/Software/miniconda3/etc/profile.d/conda.sh
conda activate cell2loc_env
export SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-32}

# ── Step 05b: Full-gene scVI integration ──────────────────────────────────────
echo "=== Step 05b: Full-gene scVI integration ==="
echo "Start: $(date)"
python "${SCRIPT_DIR}/05b_scvi_fullgene.py"
if [ $? -ne 0 ]; then echo "ERROR in Step 05b"; exit 1; fi
echo "Done:  $(date)"
echo ""

# ── Step 06: Extract urothelial cells ─────────────────────────────────────────
echo "=== Step 06: Extract urothelial cells ==="
echo "Start: $(date)"
python "${SCRIPT_DIR}/06_extract_uro.py"
if [ $? -ne 0 ]; then echo "ERROR in Step 06"; exit 1; fi
echo "Done:  $(date)"
echo ""

conda deactivate

# ── Load R environment ────────────────────────────────────────────────────────
module purge
module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0

# Prepend newer libstdc++ so zellkonverter/basilisk finds GLIBCXX_3.4.29
export LD_LIBRARY_PATH=/gpfs0/export/apps/easybuild/software/GCCcore/12.3.0/lib64:${LD_LIBRARY_PATH}

# ── Step 07: Convert h5ad → Seurat RDS ───────────────────────────────────────
echo "=== Step 07: Convert h5ad → Seurat RDS ==="
echo "Start: $(date)"
Rscript "${SCRIPT_DIR}/07_convert_to_rds.R"
if [ $? -ne 0 ]; then echo "ERROR in Step 07"; exit 1; fi
echo "Done:  $(date)"
echo ""

# ── Step 08: Harmony reclustering ─────────────────────────────────────────────
echo "=== Step 08: Harmony reclustering ==="
echo "Start: $(date)"
Rscript "${SCRIPT_DIR}/08_harmony_recluster.R"
if [ $? -ne 0 ]; then echo "ERROR in Step 08"; exit 1; fi
echo "Done:  $(date)"
echo ""

echo "======================================================"
echo "  BladderUrothelium pipeline complete: $(date)"
echo "======================================================"

sacct -j "${SLURM_JOB_ID}" \
  --format=JobID,User,Start,End,Elapsed,TotalCPU,State,NodeList \
  --noheader
