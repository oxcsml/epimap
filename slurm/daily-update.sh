#!/bin/bash
# Activate the right bash environment
source /homes/mhutchin/.bash_profile

umask 000

trap 'echo daily-update-clean: Failed before finishing with exit code $? && exit $?' ERR


# Update this repo
git pull

# Update the case data repo
cd /data/ziz/mhutchin/Rmap/covid19_datasets && git pull && cd -

# Activate the correct environment
conda activate Rmap

# Do all the data preprocessing
make preprocess-data

today=$(date +'%Y-%m-%d')
clean_directory="fits/clean-${today}"
results_directory="fits/map-${today}"
# mkdir -p $results_directory
# mkdir -p $results_directory-cori
# git rev-parse HEAD > $results_directory/git-hash.txt
# git rev-parse HEAD > $results_directory-cori/git-hash.txt

# clean 
slurm/submit-clean.sh $clean_directory

# Force recomplie to avoid mysterious bug
rm -f mapping/stan_files/Rmap.rds

# run full model
full_options="\
    --globalkernel none \
    --spatialkernel matern12 \
    --fixed_gp_time_length_scale 100.0 \
    --fixed_gp_space_length_scale 0.2 \
"
slurm/submit-run-full.sh $results_directory-full "--clean_directory $clean_directory $full_options" &
# run 2 stage model
# slurm/submit-run.sh $results_directory-two-stage "--clean_directory $clean_directory" &
# run cori model
# slurm/submit-run-cori.sh $results_directory-cori "--clean_directory $clean_directory" &
wait

# reinflate results to the website dir
dataprocessing/reinflate.sh $results_directory-full/merged_ $today &
# dataprocessing/reinflate.sh $results_directory-cori/merged_ $today-cori &
wait

# softlink to defaults
unlink docs/assets/data/default
cd docs/assets/data/ && ln -s $today default && cd -

# Update the git repo
git add docs/assets/data/$today/*
# git add docs/assets/data/$today-cori/*
git add docs/assets/data/default
git add docs/assets/data/site_data.csv
git add -f data/*
git commit -m "daily update $today"
git pull
git push


# dataprocessing/reinflate.sh $results_directory/merged_ $today &
# wait



rm -rf $results_directory/*.rds
rm -rf $clean_directory/*/*.rds
