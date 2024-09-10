#!/bin/bash
#SBATCH --array=1-128%50
#SBATCH --job-name=vpData_job
#SBATCH --partition=week-long
#SBATCH --output=vpData.out
#SBATCH --error=vpData.err

module purge
module load R

# Rscript to run an r script
# This stores which job is running (1, 2, 3, etc)
JOBID=$SLURM_ARRAY_TASK_ID
Rscript vectraPolarisdata.R $JOBID


