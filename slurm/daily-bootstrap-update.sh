#!/bin/bash


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


function cleanup()
{
  echo daily-bootstrap-update: Failed before finishing with exit code $? 
  rm -rf $results_directory/bootstrap_*/regional/*.rds
  rm -rf $results_directory/bootstrap_*/singlearea/stanfits/*.rds
  echo daily-bootstrap-update: cleaned up rds files
  exit $?
}

trap cleanup ERR

source ./slurm/cluster-config

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
slurm/submit-runs-bootstrap.sh $N_bootstrap $results_directory "${options_clean}" 1

# REGIONAL STAGE 2
slurm/submit-runs-bootstrap.sh $N_bootstrap $results_directory "${options_regional_20km}" 2

wait 

# MERGE
mkdir -p $results_directory/regional
mkdir -p $results_directory/regional/output

options_regional_20km="\
    --results_directory $results_directory \
    --approximation regional --num_regions $N_regions \
    --num_samples $N_bootstrap \
    $options_regional_20km
"
options_regional_20km="--bootstrap_merge TRUE $options_regional_20km"

echo daily-bootstrap-update: merging bootstrap samples
sbatch --wait \
    --mail-user=$USER@stats.ox.ac.uk \
    --mail-type=ALL \
    --job-name=Rmap-mergeregions_bootstrap \
    --output=$results_directory/regional/output/merge_%A_%a.out \
    --clusters=$CLUSTER \
    --partition=$PARTITION_LARGE \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=20G \
    --wrap \
    "Rscript covidmap/stage2_merge.r ${options_regional_20km}"

results_prefix="${results_directory}/regional/merged_"
dataprocessing/reinflate.sh ${results_prefix} $today

echo "copying files"
cp ${results_prefix}Rt_region.csv site_data/${today}/Rt_region.csv
cp ${results_prefix}Cpred_region.csv site_data/${today}/Cpred_region.csv
cp ${results_prefix}Cproj_region.csv site_data/${today}/Cproj_region.csv

# update website files and plots
python regional_plots/regional_plot_script.py \
            site_data/${today}/Rt_region.csv \
            site_data/${today}/Cpred_region.csv \
            site_data/${today}/Cproj_region.csv \
            docs/assets/data/region_site_data.csv \
            data/nhs_regions.csv \
            site_data/${today}

python dataprocessing/process_site_data.py

# softlink to defaults
unlink site_data/default
unlink site_data/default-bootstrap
cd site_data && ln -s $today default && cd -
cd site_data && ln -s $today default-bootstrap && cd -
rm -f docs/assets/data/default.js
echo "let map_default = \"$today\"" > docs/assets/data/default.js


# Update the git repo
git add -f site_data/${today}/*
# git add site_data/$today-cori/*
git add site_data/default
git add docs/assets/data/site_data.csv
git add docs/assets/data/region_site_data.csv
git add docs/assets/data/default.js
# git add -f data/uk_cases.csv
git commit -m "daily bootstrap update $today"
git pull
git push

rm -rf $results_directory/bootstrap_*/regional/*.rds
rm -rf $results_directory/bootstrap_*/singlearea/stanfits/*.rds

