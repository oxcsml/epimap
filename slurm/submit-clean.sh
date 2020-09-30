#!/bin/bash

set -e

if [ $# -ne 1 ]; then
  echo Usage: submit-run clean_directory
  exit 1
fi

clean_directory=$1
echo "clean_directory = $clean_directory"

options="--clean_directory $clean_directory"

mkdir $clean_directory
mkdir $clean_directory/pdfs
mkdir $clean_directory/stanfits

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
    "Rscript cleaning/clean_area.r --task_id \$SLURM_ARRAY_TASK_ID $options"
wait

echo Combining results...


sbatch --wait \
    --mail-user=$USER@stats.ox.ac.uk \
    --mail-type=ALL \
    --job-name=combine_areas \
    --output=slurm/output/cleaning/clean_combine_%A_%a.out \
    --partition=ziz-medium \
    --ntasks=1 \
    --time=18:00:00 \
    --mem-per-cpu=10G \
    --cpus-per-task=1 \
    --wrap \
    "Rscript cleaning/combine_areas.r $options"
wait
