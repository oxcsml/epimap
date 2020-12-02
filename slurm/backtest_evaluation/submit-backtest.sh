#!/bin/bash

set -e

function clean_and_map { # cleaning_directory = $1, options_clean = $2, weeks_modelled = $3, results_directory = $4, options_map = $5
    sbatch --wait \
        --mail-user=$USER@stats.ox.ac.uk \
        --mail-type=ALL \
        --job-name=clean_ts_weeks_"$3" \
        --output="$1"/output/clean_%A_%a.out \
        --partition=ziz-medium \
        --ntasks=1 \
        --time=18:00:00 \
        --mem-per-cpu=5G \
        --cpus-per-task=1 \
        --array=1-348 \
        --wrap \
        "Rscript cleaning/clean_area.r --task_id \$SLURM_ARRAY_TASK_ID $2 && echo clean_area: DONE"
    wait # are these waits needed?

    echo "submit-clean: Combining areas ($3 weeks)"

    sbatch --wait \
        --mail-user=$USER@stats.ox.ac.uk \
        --mail-type=ALL \
        --job-name=combine_areas_weeks_"$3" \
        --output="$1"/output/combine_%A_%a.out \
        --partition=ziz-medium \
        --ntasks=1 \
        --time=18:00:00 \
        --mem-per-cpu=10G \
        --cpus-per-task=1 \
        --wrap \
        "Rscript cleaning/combine_areas.r $2 && echo combine_areas: DONE"
    wait

    echo "submit-clean: DONE ($3 weeks)"

    sbatch --wait \
        --mail-user=$USER@stats.ox.ac.uk \
        --mail-type=ALL \
        --job-name=Rmap_run_weeks_"$3" \
        --output="$4"/output/run_%A_%a.out \
        --partition=ziz-large \
        --ntasks=1 \
        --cpus-per-task=1 \
        --mem-per-cpu=20G \
        --array=1-10 \
        --wrap \
        "Rscript mapping/run.r $5 \
            --cleaned_sample_id \$SLURM_ARRAY_TASK_ID \
            && echo run: DONE"

    echo submit-run: Merging results

    sbatch --wait \
        --mail-user=$USER@stats.ox.ac.uk \
        --mail-type=ALL \
        --job-name=Rmap_merge_weeks_"$3" \
        --output="$4"/output/merge_%A_%a.out \
        --partition=ziz-large \
        --ntasks=1 \
        --cpus-per-task=1 \
        --mem-per-cpu=100G \
        --wrap \
        "Rscript mapping/merge_results.r $5 && echo merge: DONE"

    wait
    
    echo submit-run: DONE

}

# TODO: turn these into positional args later?
first_day_modelled="2020-06-01"
weeks_modelled_array=(4 6 8 10 12 14 16)
backtest_directory="fits/backtests"

echo compiling
Rscript mapping/compile.r

for w in "${weeks_modelled_array[@]}"
do
    clean_directory="${backtest_directory}/clean_start_${first_day_modelled}_weeks_${w}"
    results_directory="${backtest_directory}/map_start_${first_day_modelled}_weeks_${w}"
    mkdir -p "$clean_directory"
    mkdir -p "$clean_directory"/pdfs
    mkdir -p "$clean_directory"/stanfits
    mkdir -p "$clean_directory"/output

    mkdir -p "$results_directory"/output

    options_clean="\
        --clean_directory $clean_directory \
        --first_day_modelled $first_day_modelled \
        --weeks_modelled $w \
    "

    # --time_steps 15 \ no longer needed?

    options_map="\
        --iterations 6000 \
        --observation_data cleaned_recon_sample \
        --observation_model negative_binomial_3 \
        --clean_directory $clean_directory \
        --results_directory $results_directory \
        --first_day_modelled $first_day_modelled \
        --weeks_modelled $w \
    "

    clean_and_map "$clean_directory" "$options_clean" "$w" "$results_directory" "$options_map" &

done

wait
