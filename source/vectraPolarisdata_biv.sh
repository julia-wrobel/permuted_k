#!/bin/bash
#SBATCH --array=1-128%30
#SBATCH --job-name=vpData_job
#SBATCH --mem=50G
#SBATCH --partition=wrobel,encore
#SBATCH --output=vpData.out
#SBATCH --error=vpData.err

module purge
module load R

# Rscript to run an r script
# This stores which job is running (1, 2, 3, etc)
JOBID=$SLURM_ARRAY_TASK_ID
Rscript vectraPolarisdata_biv.R $JOBID


