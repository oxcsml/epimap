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
#SBATCH --mem-per-cpu=15G # `seff' indicated that we need between 1.5-2G mem-per-cpu. 
#SBATCH --cpus-per-task=7
#
#SBATCH --array=1-6

set -e

$(sed -n "${SLURM_ARRAY_TASK_ID}p" < slurm/commands.txt)

echo 'Run completed.'

