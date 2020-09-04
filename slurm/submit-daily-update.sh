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
#SBATCH --mem-per-cpu=15G # `seff' indicated that we need between 1.5-2G mem-per-cpu. 
#SBATCH --cpus-per-task=7
#
#SBATCH --array=1

set -e

git pull

conda activap Rmap

make preprocess-data

output_file="Rmap-time-vary-reduce-latest-matern12-local-global-radiation2_uniform_in-negative_binomial_3"
time_steps=15

Rscript run-time-vary-reduce.r --time_steps=$time_steps --iterations 8000 --chains 6 --daily_update

python3 reinflate.py \
    fits/latest_updates/${output_file}-${time_steps}_Rt.csv \
    fits/latest_updates/${output_file}-${time_steps}_Pexceed.csv \
    fits/latest_updates/${output_file}-${time_steps}_Cpred.csv \
    fits/latest_updates/${output_file}-${time_steps}_Cproj.csv \
    docs/assets/data/generated/Rt.csv \
    docs/assets/data/generated/Pexceed.csv \
    docs/assets/data/generated/Cpred.csv \
    docs/assets/data/generated/Cproj.csv

git add docs/assets/data/generated/*
git commit -m "daily update"
git push


echo 'Run completed.'