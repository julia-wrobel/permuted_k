#!/bin/bash
#SBATCH --array=1-480%10
#SBATCH --job-name=var_bivariate
#SBATCH --partition=encore
#SBATCH --output=bivariate_variance.out
#SBATCH --error=bivariate_variance.err

module purge
module load R

# Rscript to run an r script
# This stores which job is running (1, 2, 3, etc)
JOBID=$SLURM_ARRAY_TASK_ID
Rscript k_bivariate_variance.R $JOBID


