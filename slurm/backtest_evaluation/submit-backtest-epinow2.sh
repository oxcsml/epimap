#!/bin/bash

set -e

cd /data/ziz/not-backed-up/scratch/szaidi/Rmap
source venv_py/bin/activate
PYTHONPATH=.

function epinow2 { # $1 = options, $2 = first_day_modelled, $3 = log_folder
    sbatch --wait \
        --mail-user=$USER@stats.ox.ac.uk \
        --mail-type=ALL \
        --job-name=Rmap_epinow2_start_"$2" \
        --output="$3"/run_%A_%a.o \
        --error="$3"/run_%A_%a.e \
        --partition=ziz-medium \
        --ntasks=1 \
        --time=18:00:00 \
        --mem-per-cpu=1G \
        --cpus-per-task="${ncores}" \
        --array=1-350%"${max_parallel_jobs}" \
        --wrap \
        "Rscript alternate_methods/epinow2_run.r --area_index \$SLURM_ARRAY_TASK_ID $1 "
    wait # are these waits needed?

    echo "Completed epinow2 over areas."

}


backtest_directory="fits/backtests_6_apr_epinow2" # relative to Rmap working directory!
first_day_modelled_array=("2020-06-29" "2020-08-10" "2020-09-07" "2020-10-05") 
weeks_modelled=15

code="/data/ziz/not-backed-up/scratch/szaidi/Rmap"
preprocess_pth=${code}/${backtest_directory}/epinow2/inputs.csv
outputs_folder=${code}/${backtest_directory}/epinow2/outputs
region_codes=${code}/${backtest_directory}/epinow2/region_codes.json
weeks_modelled=15
forecast_days=21
ncores=4 # number of cores for stan to use 
max_parallel_jobs=25

mkdir -p "${code}/${backtest_directory}/epinow2"

python alternate_methods/epinow2_preprocess.py \
    ${code}/data/uk_cases.csv \
    ${preprocess_pth} \
    --region_codes=${region_codes}

for first_day_modelled in "${first_day_modelled_array[@]}"
do
    folder="${outputs_folder}/${first_day_modelled}"
    mkdir -p $folder

    log_folder="${outputs_folder}/${first_day_modelled}/logs"
    mkdir -p $log_folder

    options="\
        --first_day_modelled ${first_day_modelled} \
        --weeks_modelled ${weeks_modelled} \
        --case_counts ${preprocess_pth} \
        --ncores ${ncores} \
        --forecast_horizon ${forecast_days} \
        --output_folder ${folder} \
        --region_codes ${region_codes} \
    "

    epinow2 "$options" "$first_day_modelled" "$log_folder" 

    python alternate_methods/epinow2_postprocess.py \
        ${outputs_folder}/${first_day_modelled} \
        ${outputs_folder} \
        ${first_day_modelled} \
        ${forecast_days} \
        --region_codes=${region_codes} \
        --prefix=${first_day_modelled}

    echo "Completed start date ${first_day_modelled}."

done

