#!/bin/bash

set -e

trap 'echo submit-run-epiestim: Failed before finishing with exit code $? && exit $?' ERR

if [ $# == 1 ]
then
  results_directory=$1
  options="\
    --results_directory $results_directory \
  "
  N=348
elif [ $# == 2 ]
then
  results_directory=$1
  options="\
    --results_directory $results_directory \
    $2
  "
  N=348
elif [ $# == 3 ]
then
  results_directory=$1
  options="\
    --results_directory $results_directory \
    $2
  "
  N=$3
else
  echo Usage: submit-run-epiestim results_directory [options] [N]
  exit 1
fi

echo "results_directory = $results_directory"

mkdir -p $results_directory
mkdir -p $results_directory/epiestim
mkdir -p $results_directory/epiestim/fits
mkdir -p $results_directory/epiestim/output

# echo submit-run-epiestim: compiling
# Rscript epimap/compile.r

echo submit-run-epiestim: running areas
sbatch --wait \
    --mail-user=$USER@stats.ox.ac.uk \
    --mail-type=ALL \
    --job-name=Rmap-epiestim \
    --output=$results_directory/epiestim/output/run_%A_%a.out \
    --partition=ziz-small \
    --ntasks=1 \
    --time=00:30:00 \
    --mem-per-cpu=5G \
    --cpus-per-task=1 \
    --array=1-$N \
    --wrap \
    "Rscript alternate_methods/epiestim_run.r --area_index \$SLURM_ARRAY_TASK_ID $options"

echo submit-run-epiestim: combining areas
sbatch --wait \
    --mail-user=$USER@stats.ox.ac.uk \
    --mail-type=ALL \
    --job-name=Rmap-combineareas \
    --output=$results_directory/epiestim/output/combine_%A_%a.out \
    --partition=ziz-small \
    --ntasks=1 \
    --time=00:30:00 \
    --mem-per-cpu=10G \
    --cpus-per-task=1 \
    --wrap \
    "Rscript alternate_methods/epiestim_combine.r $options"

wait

echo submit-run-epiestim: DONE
