#!/bin/bash
#SBATCH --array=1-120%50
#SBATCH --job-name=bothClust_variance_job
#SBATCH --partition=wrobel,encore
#SBATCH --output=bothClust_variance.out
#SBATCH --error=bothClust_variance.err

module purge
module load R

# Rscript to run an r script
# This stores which job is running (1, 2, 3, etc)
JOBID=$SLURM_ARRAY_TASK_ID
Rscript k_bothClust_variance.R $JOBID
