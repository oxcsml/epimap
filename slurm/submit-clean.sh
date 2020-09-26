#!/bin/bash

# NOTE: can replace the --dependency option with --wait for the first sbatch run if 
# this needs to be placed inside another daily running script.

areas=$(sbatch --mail-user=sheheryar.zaidi@stats.ox.ac.uk \
    --mail-type=ALL \
    --job-name=clean_ts \
    --output=slurm/output/cleaning/clean_timeseries_%A_%a.out \
    --partition=ziz-medium \
    --ntasks=1 \
    --time=18:00:00 \
    --mem-per-cpu=2G \
    --cpus-per-task=1 \
    --array=1-348 \
    --wrap \
    'Rscript preprocessing/clean_area.r --task_id $SLURM_ARRAY_TASK_ID')

prefix='Submitted batch job '
areas=${areas#"$prefix"}

sbatch --dependency=afterok:$areas \
    --mail-user=sheheryar.zaidi@stats.ox.ac.uk \
    --mail-type=ALL \
    --job-name=combine_areas \
    --output=slurm/output/cleaning/combine_areas_%A_%a.out \
    --partition=ziz-medium \
    --ntasks=1 \
    --time=18:00:00 \
    --mem-per-cpu=10G \
    --cpus-per-task=1 \
    --wrap \
    'Rscript preprocessing/combine_areas.r'

