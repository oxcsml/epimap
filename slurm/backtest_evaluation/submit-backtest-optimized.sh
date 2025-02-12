#!/bin/bash

set -e

cd /data/ziz/not-backed-up/scratch/szaidi/Rmap

function clean { # options_clean = $1, first_day_modelled = $2, results_directory = $3, options_map = $4
    sbatch --wait \
        --mail-user=$USER@stats.ox.ac.uk \
        --mail-type=ALL \
        --job-name=Rmap_singlearea_start_"$2" \
        --output="$3"/singlearea/output/run_%A_%a.out \
        --partition=ziz-medium \
        --ntasks=1 \
        --time=18:00:00 \
        --mem-per-cpu=5G \
        --cpus-per-task=1 \
        --array=1-348 \
        --wrap \
        "Rscript covidmap/stage1_run.r --area_index \$SLURM_ARRAY_TASK_ID $1 && echo clean_area: DONE"
    wait # are these waits needed?

    echo "submit-clean: Combining areas (start: $2)"

    sbatch --wait \
        --mail-user=$USER@stats.ox.ac.uk \
        --mail-type=ALL \
        --job-name=Rmap_combineareas_start_"$2" \
        --output="$3"/singlearea/output/combine_%A_%a.out \
        --partition=ziz-medium \
        --ntasks=1 \
        --time=18:00:00 \
        --mem-per-cpu=10G \
        --cpus-per-task=1 \
        --wrap \
        "Rscript covidmap/stage1_combine.r $1 && echo combine_areas: DONE"
    wait

    echo "submit-clean: DONE (start: $2)"

}

function merge { # options_clean = $1, first_day_modelled = $2, results_directory = $3, options_map = $4
    sbatch --wait \
        --mail-user=$USER@stats.ox.ac.uk \
        --mail-type=ALL \
        --job-name=Rmap_mergeregions_start_"$2" \
        --output="$3"/regional/output/merge_%A_%a.out \
        --partition=ziz-large \
        --ntasks=1 \
        --cpus-per-task=1 \
        --mem-per-cpu=10G \
        --wrap \
        "Rscript covidmap/stage2_merge.r $4 && echo merge: DONE"

    wait
    
    echo "merge: DONE"

}

# TODO: turn these into positional args for this script later?
backtest_directory="fits/backtests_5_apr_2021" # relative to Rmap working directory!

first_day_modelled_array=("2020-06-29" "2020-08-10" "2020-09-07" "2020-10-05") 
weeks_modelled=15

space_scales=("0" "0.01" "0.05" "0.1" "0.2" "0.5" "1.0")

# echo compiling
# Rscript epimap/compile.r

# create directories for all runs 
for space_scale in "${space_scales[@]}"
do
    # space_scale=${space_scales[$model_idx]}
    for first_day_modelled in "${first_day_modelled_array[@]}"
    do
        # NOTE: the signature of the results directory below 
        # should match the one in the python script used later!
        results_directory="${backtest_directory}/space_${space_scale}/start_${first_day_modelled}_weeks_${weeks_modelled}"
        mkdir -p "$results_directory"
        mkdir -p "$results_directory"/singlearea
        mkdir -p "$results_directory"/singlearea/pdfs
        mkdir -p "$results_directory"/singlearea/stanfits
        mkdir -p "$results_directory"/singlearea/output

        mkdir -p "$results_directory"/regional
        mkdir -p "$results_directory"/regional/output

        cp data/cases.csv $results_directory
        
    done

done
echo "Created all directories."

# run cleaning for all `first_day_modelled` *in parallel* (note stage1 doesn't depend on `space_scale`)
space_scale_stage1=${space_scales[0]}
for first_day_modelled in "${first_day_modelled_array[@]}"
do
    results_directory="${backtest_directory}/space_${space_scale_stage1}/start_${first_day_modelled}_weeks_${weeks_modelled}"

    options_clean="\
        --results_directory $results_directory \
        --produce_plots FALSE \
        --first_day_modelled $first_day_modelled \
        --weeks_modelled $weeks_modelled \
    "

    clean "$options_clean" "$first_day_modelled" "$results_directory" "None" &
done

wait 

# get rid of all space consuming rds files.
rm -rf /data/ziz/not-backed-up/scratch/szaidi/Rmap/${backtest_directory}/space_${space_scale_stage1}/*/singlearea/stanfits/*.rds

echo "Completed stage1 runs on slurm."

# now copy-paste the stage1 results folders for each first_day_modelled into all config folders
for space_scale in "${space_scales[@]:1}"
do 
    for first_day_modelled in "${first_day_modelled_array[@]}"
    do

    results_directory="${backtest_directory}/space_${space_scale}/start_${first_day_modelled}_weeks_${weeks_modelled}/"
    source_stage1_directory="${backtest_directory}/space_${space_scale_stage1}/start_${first_day_modelled}_weeks_${weeks_modelled}/singlearea"

    cp -R $source_stage1_directory $results_directory

    done
done

echo "Copy-pasted results from stage1 into all directories."
echo "Completed all of stage1."

# Next, use python to submit a slurm job for running the main jobs.
space_scales_list=${space_scales[0]}
for i in "${space_scales[@]:1}"; do
   space_scales_list+=,$i
done

first_day_modelled_array_list=${first_day_modelled_array[0]}
for i in "${first_day_modelled_array[@]:1}"; do
   first_day_modelled_array_list+=,$i
done

source venv_py/bin/activate
PYTHONPATH=. python slurm/backtest_evaluation/launch_backtest.py --space_scale_list=$space_scales_list --first_day_modelled_list=$first_day_modelled_array_list --backtest_directory=$backtest_directory

echo "Runs over regions (using python) completed."

# Finally, run merging like cleaning.
for space_scale in "${space_scales[@]}"
do
    for first_day_modelled in "${first_day_modelled_array[@]}"
    do
        results_directory="${backtest_directory}/space_${space_scale}/start_${first_day_modelled}_weeks_${weeks_modelled}"       
        # options_clean="\
        #     --results_directory $results_directory \
        #     --produce_plots FALSE \
        #     --first_day_modelled $first_day_modelled \
        #     --weeks_modelled $weeks_modelled \
        # "

        if [ $space_scale == 0 ]
        then
            options_map="\
                --approximation regional \
                --globalkernel none \
                --spatialkernel none \
                --results_directory $results_directory \
                --fixed_gp_time_length_scale 100.0 \
                --first_day_modelled $first_day_modelled \
                --weeks_modelled $weeks_modelled \
                --num_regions 9 \
            "
        else
            options_map="\
                --approximation regional \
                --globalkernel none \
                --spatialkernel matern12 \
                --results_directory $results_directory \
                --fixed_gp_space_length_scale $space_scale \
                --fixed_gp_time_length_scale 100.0 \
                --first_day_modelled $first_day_modelled \
                --weeks_modelled $weeks_modelled \
                --num_regions 9 \
            "
        fi

        merge "None" "$weeks_modelled" "$results_directory" "$options_map" &

    done

done

wait 

# get rid of all space consuming rds files.
rm -rf /data/ziz/not-backed-up/scratch/szaidi/Rmap/${backtest_directory}/*/*/regional/*.rds

echo "Completed."


