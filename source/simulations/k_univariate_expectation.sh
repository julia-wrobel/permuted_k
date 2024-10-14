#!/bin/bash
#SBATCH --array=1-48%48
#SBATCH --job-name=univariate_expectation_job
#SBATCH --mem=20G
#SBATCH --partition=encore
#SBATCH --output=univariate_expectation.out
#SBATCH --error=univariate_expectation.err

module purge
module load R

# Rscript to run an r script
# This stores which job is running (1, 2, 3, etc)
JOBID=$SLURM_ARRAY_TASK_ID
Rscript k_univariate_expectation.R $JOBID


