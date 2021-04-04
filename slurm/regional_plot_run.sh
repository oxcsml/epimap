#!/bin/bash

trap 'echo regional_plot_run: Failed before finishing with exit code $? && exit $?' ERR

# Activate the right bash environment
source /homes/$USER/.bash_profile
conda activate Rmap
umask 007

today=$(date +'%Y-%m-%d')
#clean_directory="fits/clean-${today}"
results_directory="fits/${today}"
mkdir -p $results_directory
cp data/cases.csv $results_directory

# clean 
slurm/submit-run-singlearea.sh $results_directory

# run full model
options_regional_20km="\
    --globalkernel none \
    --spatialkernel matern12 \
    --fixed_gp_time_length_scale 100.0 \
    --fixed_gp_space_length_scale 0.2 \
"
slurm/submit-run-regional.sh $results_directory $options_regional_20km &

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
echo completed
################### changed code ###########################
# # reinflate results to the website dir
results_prefix="${results_directory}/regional/merged_"
dataprocessing/reinflate.sh ${results_prefix} $today

cp ${results_prefix}Rt_region.csv docs/assets/data/${today}/Rt_region.csv
cp ${results_prefix}Cpred_region.csv docs/assets/data/${today}/Cpred_region.csv
cp ${results_prefix}Cproj_region.csv docs/assets/data/${today}/Cproj_region.csv

# create regional plot
python regional_plots/regional_plot_script.py \
            docs/assets/data/${today}/Rt_region.csv \
            docs/assets/data/${today}/Cpred_region.csv \
            docs/assets/data/${today}/Cproj_region.csv \
            docs/assets/data/region_site_data.csv \
            data/nhs_regions.csv \
            docs/assets/data/${today}

###########################################################
# softlink to defaults
unlink docs/assets/data/default
cd docs/assets/data/ && ln -s $today default && cd -

python dataprocessing/process_site_data.py

rm -rf $results_directory/regional/*.rds
rm -rf $results_directory/singlearea/stanfits/*.rds