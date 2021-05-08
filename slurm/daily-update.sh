#!/bin/bash

trap 'echo daily-update: Failed before finishing with exit code $? && exit $?' ERR

# Activate the right bash environment
source /homes/$USER/.profile
conda activate Rmap-daily-update
umask 007

cd /data/ziz/software/Rmap/Rmap-daily-update

# Update this repo
git pull

# Update the case data repo
cd /data/ziz/software/Rmap/covid19_datasets && git pull && cd -

python dataprocessing/process_uk_cases.py


# Activate the correct environment
# conda activate Rmap-daily-update

# Do all the data preprocessing
# make preprocess-data

today=$(date +'%Y-%m-%d')
#clean_directory="fits/clean-${today}"
results_directory="fits/${today}"
mkdir -p $results_directory
# mkdir -p $results_directory-cori
git rev-parse HEAD > $results_directory/git-hash.txt
# git rev-parse HEAD > $results_directory-cori/git-hash.txt

cp data/cases.csv $results_directory

# clean 
options_clean="\
    --weeks_modelled 20 \
    --days_ignored 7 \
"
slurm/submit-run-singlearea.sh $results_directory "$options_clean"

# Force recomplie to avoid mysterious bug
# rm -f mapping/stan_files/Rmap.rds

# run full model
options_regional_20km="\
    --globalkernel none \
    --spatialkernel matern12 \
    --fixed_gp_time_length_scale 100.0 \
    --fixed_gp_space_length_scale 0.2 \
    --weeks_modelled 15 \
    --days_ignored 7 \
    --days_predicted 2 \
    --num_steps_forcasted 3 \
"
slurm/submit-run-regional.sh $results_directory "$options_regional_20km" &

# options_regional_10km="\
#     --globalkernel none \
#     --spatialkernel matern12 \
#     --fixed_gp_time_length_scale 100.0 \
#     --fixed_gp_space_length_scale 0.1 \
# "
# slurm/submit-run-full.sh $results_directory-full-10km --clean_directory $clean_directory $options_regional_10km &

wait
# run 2 stage model
# slurm/submit-run.sh $results_directory-two-stage "--clean_directory $clean_directory" &
# run cori model
# slurm/submit-run-cori.sh $results_directory-cori "--clean_directory $clean_directory" &


# reinflate results to the website dir
dataprocessing/reinflate.sh $results_directory/regional/merged_ $today
# dataprocessing/reinflate.sh $results_directory-cori/merged_ $today-cori &

# softlink to defaults
unlink docs/assets/data/default
cd docs/assets/data/ && ln -s $today default && cd -

python dataprocessing/process_site_data.py

# Update the git repo
git add docs/assets/data/$today/*
# git add docs/assets/data/$today-cori/*
git add docs/assets/data/default
git add docs/assets/data/site_data.csv
# git add -f data/uk_cases.csv
git commit -m "daily update $today"
git pull
git push


rm -rf $results_directory/regional/*.rds
rm -rf $results_directory/singlearea/stanfits/*.rds
# rm -rf $clean_directory/*/*.rds
