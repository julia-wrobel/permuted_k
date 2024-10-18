#!/bin/bash
#SBATCH --array=1-48%20
#SBATCH --job-name=bivariate_expectation_job
#SBATCH --partition=encore
#SBATCH --output=bivariate_expectation.out
#SBATCH --error=bivariate_expectation.err

module purge
module load R

# Rscript to run an r script
# This stores which job is running (1, 2, 3, etc)
JOBID=$SLURM_ARRAY_TASK_ID
Rscript k_bivariate_expectation.R $JOBID


