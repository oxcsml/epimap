#!/bin/bash

trap 'echo daily-update: Failed before finishing with exit code $? && exit $?' ERR


CONDAROOT=/data/ziz/not-backed-up/teh/miniconda3
CONDAENVNAME=Rmap-daily-update
DIRECTORY=/data/ziz/software/Rmap/Rmap-daily-update

# Activate the right bash and conda environment
source /homes/$USER/.profile
$CONDAROOT/bin/activate 
conda activate $CONDAENVNAME
umask 007
cd $DIRECTORY

# Update this repo
git pull

# Update the case data repo
cd /data/ziz/software/Rmap/covid19_datasets && git pull && cd -

python dataprocessing/process_uk_cases.py

# Do all the data preprocessing
# make preprocess-data

if [ $# == 1 ]
then
  today=$(date +'%Y-%m-%d')-$1
else
  today=$(date +'%Y-%m-%d')
fi

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

# run full model
options_regional_10km="\
    --globalkernel none \
    --spatialkernel matern12 \
    --fixed_gp_time_length_scale 200.0 \
    --fixed_gp_space_length_scale 0.1 \
    --weeks_modelled 15 \
    --days_ignored 7 \
    --days_predicted 2 \
    --num_steps_forecasted 3 \
"
slurm/submit-run-regional.sh $results_directory "$options_regional_10km" &

# options_regional_10km="\
#     --globalkernel none \
#     --spatialkernel matern12 \
#     --fixed_gp_time_length_scale 200.0 \
#     --fixed_gp_space_length_scale 0.1 \
# "
# slurm/submit-run-full.sh $results_directory-full-10km --clean_directory $clean_directory $options_regional_10km &

wait
# run 2 stage model
# slurm/submit-run.sh $results_directory-two-stage "--clean_directory $clean_directory" &
# run cori model
# slurm/submit-run-cori.sh $results_directory-cori "--clean_directory $clean_directory" &


# reinflate results to the website dir
dataprocessing/reinflate.sh $results_directory/singlearea/ $today-singlearea
dataprocessing/reinflate.sh $results_directory/regional/merged_ $today
# dataprocessing/reinflate.sh $results_directory-cori/merged_ $today-cori &

# softlink to defaults
unlink docs/assets/data/default
cd docs/assets/data/ && ln -s $today default && cd -

python dataprocessing/process_site_data.py

# Update the git repo
git add docs/assets/data/$today/*
git add docs/assets/data/$today-singlearea/*
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
