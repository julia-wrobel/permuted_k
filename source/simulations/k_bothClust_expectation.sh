#!/bin/bash
#SBATCH --array=1-24
#SBATCH --job-name=bothClust_expectation_job
#SBATCH --partition=wrobel,encore
#SBATCH --output=bothClust_expectation.out
#SBATCH --error=bothClust_expectation.err

module purge
module load R

# Rscript to run an r script
# This stores which job is running (1, 2, 3, etc)
JOBID=$SLURM_ARRAY_TASK_ID
Rscript k_bothClust_expectation.R $JOBID
