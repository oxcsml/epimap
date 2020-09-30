#!/bin/bash

trap 'echo submit-clean: Failed before finishing with exit code $? && exit $?' ERR

echo submit-clean: Cleaning areas

sbatch --wait \
    --mail-user=$USER@stats.ox.ac.uk \
    --mail-type=ALL \
    --job-name=clean_ts \
    --output=slurm/output/cleaning/clean_%A_%a.out \
    --partition=ziz-medium \
    --ntasks=1 \
    --time=18:00:00 \
    --mem-per-cpu=5G \
    --cpus-per-task=1 \
    --array=1-348 \
    --wrap \
    'Rscript cleaning/clean_area.r \
        --task_id $SLURM_ARRAY_TASK_ID \
        && echo clean_area: DONE'

echo submit-clean: Combining areas

sbatch --wait \
    --mail-user=$USER@stats.ox.ac.uk \
    --mail-type=ALL \
    --job-name=combine_areas \
    --output=slurm/output/cleaning/combine_%A_%a.out \
    --partition=ziz-medium \
    --ntasks=1 \
    --time=18:00:00 \
    --mem-per-cpu=10G \
    --cpus-per-task=1 \
    --wrap \
    'Rscript cleaning/combine_areas.r && echo combine_areas: DONE'

echo submit-clean: DONE
