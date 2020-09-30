#!/bin/bash
# Activate the right bash environment
source /homes/mhutchin/.bash_profile

trap 'echo daily-update-clean: Failed before finishing with exit code $? && exit $?' ERR

# Update this repo
git pull

# Update the case data repo
cd /data/ziz/mhutchin/Rmap/covid19_datasets && git pull && cd -

# Activate the correct environment
conda activate Rmap

# Do all the data preprocessing
make preprocess-data
# cp -u website/site_data.csv docs/assets/data\

today=$(date +'%Y-%m-%d')
clean_directory="fits/clean-${today}"
results_directory="fits/map-${today}"
mkdir $results_directory
git rev-parse HEAD > ${results_directory}/git-hash.txt

# clean 
slurm/submit-clean.sh $clean_directory

slurm/submit-run.sh $results_directory $clean_directory

dataprocessing/reinflate.sh $results_directory $today

# soft link to latest results
rm docs/assets/data/default
ln -s docs/assets/data/$today docs/assets/data/default

# Update the git repo
git add docs/assets/data/$today/*
git add docs/assets/data/default
git add data/*
git commit -m "daily update $today"
git push
