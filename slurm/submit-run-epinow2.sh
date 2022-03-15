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
  # N=10
elif [ $# == 2 ]
then
  results_directory=$1
  options="\
    --results_directory $results_directory \
    $2
  "
  N=348
  # N=10
elif [ $# == 3 ]
then
  results_directory=$1
  options="\
    --results_directory $results_directory \
    $2
  "
  N=$3
else
  echo Usage: submit-run-epinow2 results_directory [options] [N]
  exit 1
fi

echo "options = $options"
echo "Num regions = $N"

mkdir -p $results_directory
mkdir -p $results_directory/epinow2
mkdir -p $results_directory/epinow2/stanfits
mkdir -p $results_directory/epinow2/samples
mkdir -p $results_directory/epinow2/output

echo submit-run-epinow2: running areas
sbatch --wait \
    --mail-user=$USER@stats.ox.ac.uk \
    --mail-type=ALL \
    --job-name=Rmap-epinow2 \
    --output=$results_directory/epinow2/output/run_%A_%a.out \
    --clusters=$CLUSTER \
    --partition=$PARTITION \
    --ntasks=1 \
    --time=18:00:00 \
    --mem-per-cpu=3G \
    --cpus-per-task=4 \
    --array=1-$N \
    --wrap \
    "Rscript --verbose alternate_methods/epinow2_2_run.r --area_index \$SLURM_ARRAY_TASK_ID $options"

echo submit-run-epinow2: combining areas
sbatch --wait \
    --mail-user=$USER@stats.ox.ac.uk \
    --mail-type=ALL \
    --job-name=Rmap-combineareas \
    --output=$results_directory/epinow2/output/combine_%A_%a.out \
    --clusters=$CLUSTER \
    --partition=$PARTITION \
    --ntasks=1 \
    --time=18:00:00 \
    --mem-per-cpu=10G \
    --cpus-per-task=1 \
    --wrap \
    "Rscript alternate_methods/epinow2_2_combine.r $options"

# wait

echo submit-run-epinow2: DONE
