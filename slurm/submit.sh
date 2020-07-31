#!/bin/bash
#
#SBATCH --mail-user=sheheryar.zaidi@stats.ox.ac.uk
#SBATCH --mail-type=ALL
#
#SBATCH --job-name=Rmap
#SBATCH --output=slurm/output/Rmap_ablation_%A_%a.txt
#SBATCH --partition=ziz-medium
#
#SBATCH --ntasks=1
#SBATCH --time=18:00:00
#SBATCH --mem-per-cpu=5G # `seff' indicated that we need between 1.5-2G mem-per-cpu. 
#SBATCH --cpus-per-task=6
#
#SBATCH --array=1-30

set -e

Rscript run.r --task_id $SLURM_ARRAY_TASK_ID --iterations 8000 --chains 6

echo 'Run completed.'

