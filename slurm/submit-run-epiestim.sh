#!/bin/bash

set -e

trap 'echo submit-run-epiestim: Failed before finishing with exit code $? && exit $?' ERR

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
mkdir -p $results_directory/epiestim
mkdir -p $results_directory/epiestim/pdfs
mkdir -p $results_directory/epiestim/fits
mkdir -p $results_directory/epiestim/output

cp data/cases.csv $results_directory

# echo submit-run-epiestim: compiling
# Rscript epimap/compile.r

echo submit-run-epiestim: running areas
sbatch --wait \
    --mail-user=$USER@stats.ox.ac.uk \
    --mail-type=ALL \
    --job-name=Rmap-epiestim \
    --output=$results_directory/epiestim/output/run_%A_%a.out \
    --partition=ziz-medium \
    --ntasks=1 \
    --time=18:00:00 \
    --mem-per-cpu=5G \
    --cpus-per-task=1 \
    --array=1-348 \
    --wrap \
    "Rscript alternate_methods/epiestim_run.r --area_index \$SLURM_ARRAY_TASK_ID $options && echo Rmap-epiestim: DONE"

echo submit-run-epiestim: combining areas
sbatch --wait \
    --mail-user=$USER@stats.ox.ac.uk \
    --mail-type=ALL \
    --job-name=Rmap-combineareas \
    --output=$results_directory/epiestim/output/combine_%A_%a.out \
    --partition=ziz-medium \
    --ntasks=1 \
    --time=18:00:00 \
    --mem-per-cpu=10G \
    --cpus-per-task=1 \
    --wrap \
    "Rscript alternate_methods/epiestim_combine.r $options && echo Rmap-combineareas: DONE"

wait

echo submit-run-epiestim: DONE
