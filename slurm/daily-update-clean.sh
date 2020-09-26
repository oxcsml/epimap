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

# Submit each region to clean 
slurm/submit-clean.sh

# Submit the job and wait for completion
sbatch --wait slurm/submit-daily-update-cleaned.sh 
wait

# Recombine samples
results_directory="fits/Rmap-cleaned-$(date +'%Y-%m-%d')"

Rscript postprocess_samples.r --time_steps=15 \
    --iterations 100 \
    --chains 1 \
    --results_directory $results_directory \

# Reinflate the results to the website directory
python3 reinflate.py \
    ${results_directory}/merged_Rt.csv \
    ${results_directory}/merged_Pexceed.csv \
    ${results_directory}/merged_Cweekly.csv \
    ${results_directory}/merged_Cpred.csv \
    ${results_directory}/merged_Cproj.csv \
    docs/assets/data/default/Rt.csv \
    docs/assets/data/default/Pexceed.csv \
    docs/assets/data/default/Cweekly.csv \
    docs/assets/data/default/Cpred.csv \
    docs/assets/data/default/Cproj.csv

# Update the git repo
git add docs/assets/data/default/*
git add data/*
git commit -m "daily update"
git push
