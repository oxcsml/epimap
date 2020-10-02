#!/bin/bash

set -e

trap 'echo submit-clean: Failed before finishing with exit code $? && exit $?' ERR

echo submit-clean: Cleaning areas

if [ $# -ne 1 ]; then
  echo Usage: submit-run clean_directory
  exit 1
fi

clean_directory=$1
echo "clean_directory = $clean_directory"

options="--clean_directory $clean_directory"

mkdir -p $clean_directory
mkdir -p $clean_directory/pdfs
mkdir -p $clean_directory/stanfits
mkdir -p $clean_directory/output

sbatch --wait \
    --mail-user=$USER@stats.ox.ac.uk \
    --mail-type=ALL \
    --job-name=clean_ts \
    --output=$clean_directory/output/clean_%A_%a.out \
    --partition=ziz-medium \
    --ntasks=1 \
    --time=18:00:00 \
    --mem-per-cpu=5G \
    --cpus-per-task=1 \
    --array=1-348 \
    --wrap \
    "Rscript cleaning/clean_area.r --task_id \$SLURM_ARRAY_TASK_ID $options && echo clean_area: DONE"
wait

echo submit-clean: Combining areas

sbatch --wait \
    --mail-user=$USER@stats.ox.ac.uk \
    --mail-type=ALL \
    --job-name=combine_areas \
    --output=$clean_directory/output/combine_%A_%a.out \
    --partition=ziz-medium \
    --ntasks=1 \
    --time=18:00:00 \
    --mem-per-cpu=10G \
    --cpus-per-task=1 \
    --wrap \
    "Rscript cleaning/combine_areas.r $options && echo combine_areas: DONE"
wait

echo submit-clean: DONE
