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
# cp -u website/site_data.csv docs/assets/data\

results_directory="fits/Rmap-cleaned-$(date +'%Y-%m-%d')"
mkdir $results_directory
git rev-parse HEAD > ${results_directory}/git-hash.txt

# Submit each region to clean and smooth
sbatch --wait \
    --mail-user=michael.hutchinson@stats.ox.ac.uk \
    --mail-type=ALL \
    --job-name=clean_ts \
    --output=slurm/output/cleaning/clean_timeseries_%A_%a.out \
    --partition=ziz-medium \
    --ntasks=1 \
    --time=18:00:00 \
    --mem-per-cpu=5G \
    --cpus-per-task=1 \
    --array=1-348 \
    --wrap \
    'Rscript preprocessing/clean_area.r --task_id $SLURM_ARRAY_TASK_ID'
wait

# Recombine samples from regions into country samples
sbatch --wait \
    --mail-user=michael.hutchinson@stats.ox.ac.uk \
    --mail-type=ALL \
    --job-name=combine_areas \
    --output=slurm/output/cleaning/combine_areas_%A_%a.out \
    --partition=ziz-medium \
    --ntasks=1 \
    --time=18:00:00 \
    --mem-per-cpu=10G \
    --cpus-per-task=1 \
    --wrap \
    'Rscript preprocessing/combine_areas.r'
wait

# Submit the job and wait for completion
sbatch --wait slurm/submit-daily-update-cleaned.sh 
wait

# Recombine samples

Rscript postprocess_samples.r --time_steps=15 \
    --iterations 8000 \
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
git pull
git add docs/assets/data/default/*
git add docs/assets/data/*
git add data/*
git commit -m "daily update"
git push
