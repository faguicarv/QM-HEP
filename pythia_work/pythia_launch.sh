#!/bin/bash

#SBATCH --job-name=pythia_Z-decay_tomography
#SBATCH --output=logs/pythia_%A_%a.out
#SBATCH --error=logs/pythia_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00
#SBATCH --mem=2G
#SBATCH --array=1-10

mkdir -p logs

module load pythia/8.315

./e+e- $SLURM_ARRAY_TASK_ID
