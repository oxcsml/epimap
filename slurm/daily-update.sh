#!/bin/bash
# Activate the right bash environment
source /homes/mhutchin/.bash_profile

umask 000

trap 'echo daily-update-clean: Failed before finishing with exit code $? && exit $?' ERR

rm -rf map-*/*.rds
rm -rf clean-/**/*.rds

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
mkdir -p $results_directory
# mkdir -p $results_directory-cori
git rev-parse HEAD > $results_directory/git-hash.txt
# git rev-parse HEAD > $results_directory-cori/git-hash.txt

# clean 
slurm/submit-clean.sh $clean_directory

rm -f mapping/stan_files/Rmap.rds

dataprocessing/reinflate.sh $clean_directory/ $today &
# dataprocessing/reinflate.sh $results_directory-cori/merged_ $today-cori &
wait

# softlink to defaults
unlink docs/assets/data/default
cd docs/assets/data/ && ln -s $today default && cd -

# Update the git repo
git add docs/assets/data/$today/*
# git add docs/assets/data/$today-cori/*
git add docs/assets/data/default
git add -f data/*
git commit -m "daily update $today"
git push


# dataprocessing/reinflate.sh $results_directory/merged_ $today &
# wait

slurm/submit-run.sh $results_directory $clean_directory &
# slurm/submit-run-cori.sh $results_directory-cori $clean_directory &
wait