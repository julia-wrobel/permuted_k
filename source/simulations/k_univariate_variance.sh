#!/bin/bash
#SBATCH --array=1-1600%40
#SBATCH --job-name=univariate_variance_job
#SBATCH --mem=50G
#SBATCH --partition=wrobel
#SBATCH --output=univariate_variance.out
#SBATCH --error=univariate_variance.err

module purge
module load R

# Rscript to run an r script
# This stores which job is running (1, 2, 3, etc)
JOBID=$SLURM_ARRAY_TASK_ID
Rscript k_univariate_variance.R $JOBID


