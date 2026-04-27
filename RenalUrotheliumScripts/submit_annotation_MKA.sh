#!/bin/bash
#SBATCH --job-name=RenalUro_annot_MKA
#SBATCH --output=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/RenalUrotheliumScripts/logs/annotation_MKA_%j.out
#SBATCH --error=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/RenalUrotheliumScripts/logs/annotation_MKA_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --mem=600G
#SBATCH --time=72:00:00
#SBATCH --partition=himem
# To use a GPU node instead, comment out the two lines above and uncomment:
# #SBATCH --partition=gpu
# #SBATCH --gres=gpu:1
# #SBATCH --mem=300G

SCRIPT_DIR="/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/UsedSingleCells/RenalUrotheliumScripts"

mkdir -p "${SCRIPT_DIR}/logs"
mkdir -p "${SCRIPT_DIR}/output/annotation_MKA"

echo "======================================================"
echo "  RenalUrothelium — MKA annotation"
echo "  Host  : $(hostname)"
echo "  Start : $(date)"
echo "======================================================"

source /home/gdbecknelllab/xxw004/Software/miniconda3/etc/profile.d/conda.sh
conda activate cell2loc_env

echo "Python:     $(python --version)"
echo "scvi-tools: $(python -c 'import scvi; print(scvi.__version__)')"
echo ""

python "${SCRIPT_DIR}/annotation_MKA.py"

EXIT_CODE=$?
echo ""
echo "End time: $(date)"
echo "Exit code: ${EXIT_CODE}"

sacct -j "${SLURM_JOB_ID}" \
  --format=JobID,User,Start,End,Elapsed,TotalCPU,MaxRSS,State \
  --noheader

exit ${EXIT_CODE}
