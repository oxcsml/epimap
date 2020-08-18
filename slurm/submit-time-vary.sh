#!/bin/bash
#
#SBATCH --mail-user=michael.hutchinson@stats.ox.ac.uk
#SBATCH --mail-type=ALL
#
#SBATCH --job-name=Rmap-time-vary
#SBATCH --output=slurm/output/Rmap_ablation_%A_%a.txt
#SBATCH --partition=ziz-large
#
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=5G # `seff' indicated that we need between 1.5-2G mem-per-cpu. 
#SBATCH --cpus-per-task=6
#
#SBATCH --array=1-30

set -e

Rscript run-time-vary.r --time_steps=23 --iterations 8000 --chains 6  --task_id $SLURM_ARRAY_TASK_ID

echo 'Run completed.'

