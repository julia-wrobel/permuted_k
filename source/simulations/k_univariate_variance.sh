#!/bin/bash
#SBATCH --array=1-528%10
#SBATCH --job-name=univariate_variance_job
#SBATCH --partition=wrobel,encore
#SBATCH --output=univariate_variance.out
#SBATCH --error=univariate_variance.err

module purge
module load R

# Rscript to run an r script
# This stores which job is running (1, 2, 3, etc)
JOBID=$SLURM_ARRAY_TASK_ID
Rscript k_univariate_variance.R $JOBID


