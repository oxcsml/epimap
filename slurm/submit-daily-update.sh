#!/bin/bash
#
#SBATCH --mail-user=michael.hutchinson@stats.ox.ac.uk
#SBATCH --mail-type=ALL
#
#SBATCH --job-name=Rmap-daily-update
#SBATCH --output=slurm/output/Rmap_daily-update_%A_%a.txt
#SBATCH --partition=ziz-large
#
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=20G # `seff' indicated that we need between 1.5-2G mem-per-cpu. 
#SBATCH --cpus-per-task=7
#
#SBATCH --array=1

set -e

Rscript run-time-vary-reduce.r --time_steps=15 --iterations 8000 --chains 6 --daily_update

echo 'Run completed.'