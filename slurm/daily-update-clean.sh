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

jobname=$(date +'%Y-%m-%d')

# clean 
slurm/submit-clean.sh

slurm/submit-run.sh $jobname

# soft link to latest results
rm docs/assets/data/default
ln -s docs/assets/data/$jobname docs/assets/data/default

# Update the git repo
git add docs/assets/data/$jobname/*
git add docs/assets/data/default
git add data/*
git commit -m "daily update $jobname"
git push
