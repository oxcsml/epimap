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
  today=run-all-$1
else
  today=run-all
fi

#clean_directory="fits/clean-${today}"
results_directory="fits/${today}"
mkdir -p $results_directory
# mkdir -p $results_directory-cori
# git rev-parse HEAD > $results_directory/git-hash.txt
# git rev-parse HEAD > $results_directory-cori/git-hash.txt

cp data/cases.csv $results_directory

# clean 
options_clean="\
    --weeks_modelled 25 \
    --last_day_modelled 2021-01-30  \
    --days_predicted 2 \
    --num_steps_forecasted 3 \
"
slurm/submit-run-singlearea.sh $results_directory "$options_clean"

# run full model
options_regional_20km="\
    --globalkernel none \
    --spatialkernel matern12 \
    --fixed_gp_time_length_scale 200.0 \
    --fixed_gp_space_length_scale 0.2 \
    --weeks_modelled 22 \
    --last_day_modelled 2021-01-30  \
    --days_predicted 2 \
    --num_steps_forecasted 3 \
"
slurm/submit-run-regional.sh $results_directory "$options_regional_20km" &

wait


# reinflate results to the website dir
dataprocessing/reinflate.sh $results_directory/singlearea/ $today-singlearea
dataprocessing/reinflate.sh $results_directory/regional/merged_ $today
# dataprocessing/reinflate.sh $results_directory-cori/merged_ $today-cori &

# softlink to defaults
#unlink docs/assets/data/default
#cd docs/assets/data/ && ln -s $today default && cd -

python dataprocessing/process_site_data.py

# Update the git repo
#git add docs/assets/data/$today/*
#git add docs/assets/data/$today-singlearea/*
# git add docs/assets/data/$today-cori/*
#git add docs/assets/data/default
#git add docs/assets/data/site_data.csv
# git add -f data/uk_cases.csv
#git commit -m "daily update $today"
#git pull
#git push


#rm -rf $results_directory/regional/*.rds
#rm -rf $results_directory/singlearea/stanfits/*.rds
# rm -rf $clean_directory/*/*.rds
