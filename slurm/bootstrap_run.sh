#!/bin/bash
trap 'echo regional_plot_run: Failed before finishing with exit code $? && exit $?' ERR

# Activate the right bash environment
source /homes/$USER/.bashrc
conda activate Rmap
umask 007

today=$(date +'%Y-%m-%d')
#clean_directory="fits/clean-${today}"
results_directory="fits/${today}"

options_clean="\
    --weeks_modelled 15 \
    --days_ignored 7 \
"
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
if [ $# == 0 ]
then
    options_clean="\
        --Aip 2.29 --Adp 1.57 --Bip 0.36 --Bdp 0.65 \
        ${options_clean}
    "
    options_regional_20km="\
        --Aip 2.29 --Adp 1.57 --Bip 0.36 --Bdp 0.65 \
        ${options_regional_20km}
    "
elif [ $# == 4 ]
then
    bootstrap_id=$1
    results_directory=$2
    options_clean=$3
    options_clean="\
        $(sed -n "${bootstrap_id}p" < ${results_directory}/bootstrap_params.txt) \
        ${options_clean}
    "
    options_regional_20km=$4
    options_regional_20km="\
        $(sed -n "${bootstrap_id}p" < ${results_directory}/bootstrap_params.txt) \
        ${options_regional_20km}
    "
    results_directory="${results_directory}/bootstrap_$1"
else
  echo Usage: regions_as_areas_stage1_run [bootstrap_id results_directory options_clean]
  exit 1
fi
mkdir -p $results_directory
cp data/cases.csv $results_directory

# clean 
slurm/submit-run-singlearea.sh $results_directory "${options_clean}"

# regional
slurm/submit-run-regional.sh $results_directory "${options_regional_20km}" &

wait
echo completed
################### changed code ###########################
# # reinflate results to the website dir
# results_prefix="${results_directory}/regional/merged_"
# dataprocessing/reinflate.sh ${results_prefix} $today

# cp ${results_prefix}Rt_region.csv docs/assets/data/${today}/Rt_region.csv
# cp ${results_prefix}Cpred_region.csv docs/assets/data/${today}/Cpred_region.csv
# cp ${results_prefix}Cproj_region.csv docs/assets/data/${today}/Cproj_region.csv

# # create regional plot
# python regional_plots/regional_plot_script.py \
#             docs/assets/data/${today}/Rt_region.csv \
#             docs/assets/data/${today}/Cpred_region.csv \
#             docs/assets/data/${today}/Cproj_region.csv \
#             docs/assets/data/region_site_data.csv \
#             data/nhs_regions.csv \
#             docs/assets/data/${today}/regional_plot.pdf

