#!/bin/bash
#SBATCH --job-name=RenalUro_annot
#SBATCH --output=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/RenalUrotheliumScripts/logs/annot_%j.out
#SBATCH --error=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/RenalUrotheliumScripts/logs/annot_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=6:00:00
#SBATCH --partition=himem

SCRIPT_DIR="/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/RenalUrotheliumScripts"

mkdir -p "${SCRIPT_DIR}/logs"

echo "======================================================"
echo "  RenalUrothelium: Add Lake + MKA Annotations"
echo "  Host  : $(hostname)"
echo "  Start : $(date)"
echo "======================================================"

# ── Step 7a: Add annotations to h5ad (Python) ─────────────────────────────────
source /home/gdbecknelllab/xxw004/Software/miniconda3/etc/profile.d/conda.sh
conda activate cell2loc_env

echo "=== Step 7a: Add Lake + MKA annotations → h5ad ==="
echo "Start: $(date)"
python "${SCRIPT_DIR}/07_add_annotations.py"
if [ $? -ne 0 ]; then echo "ERROR in Step 7a"; exit 1; fi
echo "Done:  $(date)"
echo ""

conda deactivate

# ── Step 7b: Convert annotated h5ad → Seurat RDS (R) ─────────────────────────
module purge
module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0

# Prepend newer libstdc++ so zellkonverter/basilisk finds GLIBCXX_3.4.29
export LD_LIBRARY_PATH=/gpfs0/export/apps/easybuild/software/GCCcore/12.3.0/lib64:${LD_LIBRARY_PATH}

echo "=== Step 7b: Convert annotated h5ad → Seurat RDS ==="
echo "Start: $(date)"
Rscript "${SCRIPT_DIR}/07b_convert_annotated_to_rds.R"
if [ $? -ne 0 ]; then echo "ERROR in Step 7b"; exit 1; fi
echo "Done:  $(date)"
echo ""

echo "======================================================"
echo "  Done: $(date)"
echo "======================================================"

sacct -j "${SLURM_JOB_ID}" \
  --format=JobID,User,Start,End,Elapsed,TotalCPU,State,NodeList \
  --noheader
