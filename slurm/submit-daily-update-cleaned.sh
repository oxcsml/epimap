#!/bin/bash
#
#SBATCH --mail-user=$USER@stats.ox.ac.uk
#SBATCH --mail-type=ALL
#
#SBATCH --job-name=Rmap-daily-update-clean
#SBATCH --output=slurm/output/Rmap_daily-update-cleaned_%A_%a.txt
#SBATCH --partition=ziz-large
#
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G # `seff' indicated that we need between 1.5-2G mem-per-cpu. 
#
#SBATCH --array=1-10

set -e

results_directory="fits/Rmap-cleaned-$(date +'%Y-%m-%d')"

Rscript run.r --time_steps=15 \
    --iterations 8000 \
    --chains 1 \
    --cleaned_sample_id ${SLURM_ARRAY_TASK_ID} \
    --results_directory $results_directory \

echo 'Run completed.'