#!/bin/bash
#SBATCH --job-name=AllUro_integrate
#SBATCH --output=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/FinalUrotheliumCells/UrotheliumScripts/logs/AllUro_integrate_%j.out
#SBATCH --error=/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/FinalUrotheliumCells/UrotheliumScripts/logs/AllUro_integrate_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=1000G
#SBATCH --time=6:00:00
#SBATCH --partition=himem

SCR_DIR="/vast0/home/gdjacksonlab/lab/xxw004/UUO/Datasets/Mouse/FinalUrotheliumCells/UrotheliumScripts"
mkdir -p "${SCR_DIR}/logs"

module purge
module load GCC/9.3.0 OpenMPI/4.0.3 R/4.4.0

echo "Hostname  : $(hostname)"
echo "Start time: $(date)"
echo ""

Rscript "${SCR_DIR}/01_integrate_AllUrothelium_harmony.R"

echo ""
echo "End time: $(date)"
sacct -j "${SLURM_JOB_ID}" \
  --format=JobID,User,Start,End,Elapsed,TotalCPU,State,NodeList \
  --noheader --units=G
