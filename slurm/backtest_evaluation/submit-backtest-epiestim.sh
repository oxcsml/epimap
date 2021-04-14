#!/bin/bash

set -e

cd /data/ziz/not-backed-up/scratch/szaidi/Rmap

function epiestim { # options = $1, first_day_modelled = $2, results_directory = $3
    echo submit-run-epiestim: running areas
    sbatch --wait \
        --mail-user=$USER@stats.ox.ac.uk \
        --mail-type=ALL \
        --job-name=Rmap_epiestim_start_"$2" \
        --output="$3"/epiestim/output/run_%A_%a.out \
        --partition=ziz-medium \
        --ntasks=1 \
        --time=18:00:00 \
        --mem-per-cpu=5G \
        --cpus-per-task=1 \
        --array=1-348 \
        --wrap \
        "Rscript alternate_methods/epiestim_run.r --area_index \$SLURM_ARRAY_TASK_ID $1"

    echo submit-run-epiestim: combining areas
    sbatch --wait \
        --mail-user=$USER@stats.ox.ac.uk \
        --mail-type=ALL \
        --job-name=Rmap_epiestim_combineareas_start_"$2" \
        --output="$3"/epiestim/output/combine_%A_%a.out \
        --partition=ziz-medium \
        --ntasks=1 \
        --time=18:00:00 \
        --mem-per-cpu=10G \
        --cpus-per-task=1 \
        --wrap \
        "Rscript alternate_methods/epiestim_combine.r $1"

    wait

    echo "submit-run-epiestim on start $2: DONE"

}


backtest_directory="fits/backtests_8_apr_2021_epiestim" # relative to Rmap working directory! 
first_day_modelled_array=("2020-06-29" "2020-08-10" "2020-09-07" "2020-10-05") 
weeks_modelled=15


for first_day_modelled in "${first_day_modelled_array[@]}"
do
    results_directory="${backtest_directory}/epiestim/start_${first_day_modelled}_weeks_${weeks_modelled}"

    mkdir -p "$results_directory"
    mkdir -p "$results_directory"/epiestim
    mkdir -p "$results_directory"/epiestim/pdfs
    mkdir -p "$results_directory"/epiestim/fits
    mkdir -p "$results_directory"/epiestim/output

    cp data/cases.csv $results_directory

    options="\
        --results_directory $results_directory \
        --produce_plots FALSE \
        --first_day_modelled $first_day_modelled \
        --weeks_modelled $weeks_modelled \
    "

    epiestim "$options" "$first_day_modelled" "$results_directory" &
    
done

wait 

