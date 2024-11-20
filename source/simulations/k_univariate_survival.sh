#!/bin/bash
#SBATCH --array=1-200%50
#SBATCH --job-name=surv_job
#SBATCH --partition=encore
#SBATCH --output=surv.out
#SBATCH --error=surv.err

module purge
module load R

# Rscript to run an r script
# This stores which job is running (1, 2, 3, etc)
JOBID=$SLURM_ARRAY_TASK_ID
Rscript k_univariate_survival.R $JOBID


