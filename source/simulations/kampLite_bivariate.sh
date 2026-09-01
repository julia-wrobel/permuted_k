#!/bin/bash
#SBATCH --array=1-54%50
#SBATCH --job-name=bivariate_kamplite_job
#SBATCH --partition=wrobel,encore
#SBATCH --output=bivariate_kamplite.out
#SBATCH --error=bivariate_kamplite.err

module purge
module load R

# Rscript to run an r script
# This stores which job is running (1, 2, 3, etc)
JOBID=$SLURM_ARRAY_TASK_ID
Rscript kampLite_bivariate.R $JOBID

