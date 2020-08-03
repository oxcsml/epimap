#!/bin/bash
#
#SBATCH --mail-user=sheheryar.zaidi@stats.ox.ac.uk
#SBATCH --mail-type=ALL
#
#SBATCH --job-name=Rmap-time-vary
#SBATCH --output=slurm/output/Rmap_ablation_%A_%a.txt
#SBATCH --partition=ziz-medium
#
#SBATCH --ntasks=1
#SBATCH --time=18:00:00
#SBATCH --mem-per-cpu=5G # `seff' indicated that we need between 1.5-2G mem-per-cpu. 
#SBATCH --cpus-per-task=6
#
#SBATCH --array=1

set -e

Rscript run-time-vary.r --time_steps=10 --iterations 8000 --chains 6 --task_id $SLURM_ARRAY_TASK_ID

mv fits/* ../../not-backed-up/mhutchin/Rmap/fits/

echo 'Run completed.'
