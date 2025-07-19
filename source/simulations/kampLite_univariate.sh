#!/bin/bash
#SBATCH --array=1-24%24
#SBATCH --job-name=univariate_kamplite_job
#SBATCH --partition=encore
#SBATCH --output=univariate_kamplite.out
#SBATCH --error=univariate_kamplite.err

module purge
module load R

# Rscript to run an r script
# This stores which job is running (1, 2, 3, etc)
JOBID=$SLURM_ARRAY_TASK_ID
Rscript kampLite_univariate.R $JOBID


