#!/bin/bash
# Activate the right bash environment
source /homes/mhutchin/.bash_profile

# Update this repo
git pull

# Update the case data repo
cd /data/ziz/mhutchin/Rmap/covid19_datasets && git pull && cd -

# Activate the correct environment
conda activate Rmap

# Do all the data preprocessing
make preprocess-data
cp -u website/site_data.csv docs/assets/data

output_file="Rmap-time-vary-reduce-latest-matern12-local-global-radiation2_uniform_in-negative_binomial_3"
time_steps=15

# Submit the job and wait for completion
sbatch --wait slurm/submit-daily-update.sh 
wait

# Reinflate the results to the website directory
python3 reinflate.py \
    fits/latest_updates/${output_file}-${time_steps}_Rt.csv \
    fits/latest_updates/${output_file}-${time_steps}_Pexceed.csv \
    fits/latest_updates/${output_file}-${time_steps}_Cweekly.csv \
    fits/latest_updates/${output_file}-${time_steps}_Cpred.csv \
    fits/latest_updates/${output_file}-${time_steps}_Cproj.csv \
    docs/assets/data/generated/Rt.csv \
    docs/assets/data/generated/Pexceed.csv \
    docs/assets/data/generated/Cweekly.csv \
    docs/assets/data/generated/Cpred.csv \
    docs/assets/data/generated/Cproj.csv

# Update the git repo
git add docs/assets/data/generated/*
git add data/*
git commit -m "daily update"
git push
