#!/bin/bash
trap 'echo daily-bootstrap-update: Failed before finishing with exit code $? && exit $?' ERR

CONDAROOT=/data/ziz/not-backed-up/teh/miniconda3
CONDAENVNAME=Rmap-daily-update
DIRECTORY=/data/ziz/software/Rmap/Rmap-daily-update

# Activate the right bash environment
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
python dataprocessing/process_region_site_data.py 

if [ $# == 1 ]
then
  today=$(date +'%Y-%m-%d')-bootstrap-$1
else
  today=$(date +'%Y-%m-%d')-bootstrap
fi

results_directory="fits/${today}"

mkdir -p $results_directory
git rev-parse HEAD > $results_directory/git-hash.txt

cp data/cases.csv $results_directory

options_clean="\
    --weeks_modelled 15 \
    --days_ignored 7 \
"

# run full model
options_regional_20km="\
    --globalkernel none \
    --spatialkernel matern12 \
    --fixed_gp_time_length_scale 100.0 \
    --fixed_gp_space_length_scale 0.2 \
    --weeks_modelled 15 \
    --days_ignored 7 \
    --days_predicted 2 \
    --num_steps_forecasted 3 \
"

N_bootstrap=10
N_regions=9

# CREATE SAMPLES
python3 regional_plots/create_bootstrap_samples.py --save_dir $results_directory --num_samples $N_bootstrap

# CLEAN STAGE 1
slurm/submit-bootstrap-runs-singlearea.sh $N_bootstrap $results_directory "${options_clean}"

# REGIONAL STAGE 2
slurm/submit-bootstrap-runs-regional.sh $N_bootstrap $results_directory "${options_regional_20km}"

wait 


results_prefix="${results_directory}/regional/merged_"
dataprocessing/reinflate.sh ${results_prefix} $today

# softlink to defaults
unlink docs/assets/data/default && unlink docs/assets/data/default-bootstrap
cd docs/assets/data/ && ln -s $today default && ln -s $today default-bootstrap && cd -

echo "Copying files"
cp ${results_prefix}Rt_region.csv docs/assets/data/${today}/Rt_region.csv
cp ${results_prefix}Cpred_region.csv docs/assets/data/${today}/Cpred_region.csv
cp ${results_prefix}Cproj_region.csv docs/assets/data/${today}/Cproj_region.csv

# update website files and plots
python regional_plots/regional_plot_script.py \
            docs/assets/data/${today}/Rt_region.csv \
            docs/assets/data/${today}/Cpred_region.csv \
            docs/assets/data/${today}/Cproj_region.csv \
            docs/assets/data/region_site_data.csv \
            data/nhs_regions.csv \
            docs/assets/data/${today}

python dataprocessing/process_site_data.py

# Update the git repo
git add -f docs/assets/data/${today}/*
git add docs/assets/data/default-bootstrap
git add docs/assets/data/default
git add docs/assets/data/site_data.csv
git add docs/assets/data/region_site_data.csv
# git add -f data/uk_cases.csv
git commit -m "weekly regional update $today"
git pull
git push

echo "Deleting Stan fits"
rm -rf $results_directory/bootstrap_*/regional/*.rds
rm -rf $results_directory/bootstrap_*/singlearea/stanfits/*.rds

