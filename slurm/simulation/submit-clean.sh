#!/bin/bash

set -e

trap 'echo submit-clean: Failed before finishing with exit code $? && exit $?' ERR

if [ $# -ne 3 ]; then
  echo Usage: submit-clean clean_directory data_directory n_regions
  exit 1
fi

clean_directory=$1
data_directory=$2
n_regions=$3
echo "clean_directory = $clean_directory"
echo "data_directory = $data_directory"

options="\
--clean_directory $clean_directory \
--data_directory $data_directory \
--weeks_modelled 25 \
"

mkdir -p $clean_directory
mkdir -p $clean_directory/pdfs
mkdir -p $clean_directory/stanfits
mkdir -p $clean_directory/output

echo submit-clean: compiling
Rscript epimap/compile.r

echo submit-clean: cleaning area
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
    --array=1-$n_regions \
    --wrap \
    "Rscript covidmap/stage1_run.r --area_index \$SLURM_ARRAY_TASK_ID $options && echo clean_area: DONE"
wait

echo submit-clean: combining areas
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
    "Rscript covidmap/stage1_combine.r $options && echo combine_areas: DONE"
wait

echo submit-clean: DONE
