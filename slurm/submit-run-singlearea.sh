#!/bin/bash

set -e

trap 'echo submit-clean: Failed before finishing with exit code $? && exit $?' ERR

if [ $# -ne 1 ]; then
  echo Usage: submit-run results_directory
  exit 1
fi

results_directory=$1
echo "results_directory = $results_directory"

options="\
  --produce_plots TRUE \
  --results_directory $results_directory \
"

mkdir -p $results_directory
mkdir -p $results_directory/singlearea
mkdir -p $results_directory/singlearea/pdfs
mkdir -p $results_directory/singlearea/stanfits
mkdir -p $results_directory/singlearea/output

cp data/cases.csv $results_directory

echo submit-run-singlearea: compiling
Rscript epimap/compile.r

echo submit-run-singlearea: running areas
sbatch --wait \
    --mail-user=$USER@stats.ox.ac.uk \
    --mail-type=ALL \
    --job-name=clean_ts \
    --output=$results_directory/singlearea/output/clean_%A_%a.out \
    --partition=ziz-medium \
    --ntasks=1 \
    --time=18:00:00 \
    --mem-per-cpu=5G \
    --cpus-per-task=1 \
    --array=1-348 \
    --wrap \
    "Rscript covidmap/stage1_run.r --area_index \$SLURM_ARRAY_TASK_ID $options && echo clean_area: DONE"

wait

echo submit-run-singlearea: combining areas
sbatch --wait \
    --mail-user=$USER@stats.ox.ac.uk \
    --mail-type=ALL \
    --job-name=combine_areas \
    --output=$results_directory/singlearea/output/combine_%A_%a.out \
    --partition=ziz-medium \
    --ntasks=1 \
    --time=18:00:00 \
    --mem-per-cpu=10G \
    --cpus-per-task=1 \
    --wrap \
    "Rscript covidmap/stage1_combine.r $options && echo combine_areas: DONE"

wait

echo submit-clean: DONE
