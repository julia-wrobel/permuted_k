#!/bin/bash
#SBATCH --array=1-1000%50
#SBATCH --job-name=var_perm_job
#SBATCH --mem=50G
#SBATCH --partition=encore
#SBATCH --output=var_perm.out
#SBATCH --error=var_perm.err

module purge
module load R

# Rscript to run an r script
# This stores which job is running (1, 2, 3, etc)
JOBID=$SLURM_ARRAY_TASK_ID
Rscript k_univariate_variance_permOnly.R $JOBID


