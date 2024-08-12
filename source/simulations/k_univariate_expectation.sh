#!/bin/bash
#SBATCH --array=1-24%24
#SBATCH --job-name=univariate_expectation_job
#SBATCH --partition=short-cpu
#SBATCH --output=univariate_expectation.out
#SBATCH --error=univariate_expectation.err

module purge
module load R

# Rscript to run an r script
# This stores which job is running (1, 2, 3, etc)
JOBID=$SLURM_ARRAY_TASK_ID
Rscript k_univariate_expectation.R $JOBID


