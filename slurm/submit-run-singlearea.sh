#!/bin/bash

set -e

trap 'echo submit-run-singlearea: Failed before finishing with exit code $? && exit $?' ERR

if [ $# == 1 ]
then
  results_directory=$1
  options="\
    --produce_plots FALSE \
    --results_directory $results_directory \
  "
  N=348
elif [ $# == 2 ]
then
  results_directory=$1
  options="\
    --produce_plots FALSE \
    --results_directory $results_directory \
    $2
  "
  N=348
elif [ $# == 3 ]
then
  results_directory=$1
  options="\
    --produce_plots FALSE \
    --results_directory $results_directory \
    $2
  "
  N=$3
else
  echo Usage: submit-run-singlearea results_directory \"[options]\" [N]
  exit 1
fi

echo "results_directory = $results_directory"
echo $options

mkdir -p $results_directory
mkdir -p $results_directory/singlearea
mkdir -p $results_directory/singlearea/pdfs
mkdir -p $results_directory/singlearea/stanfits
mkdir -p $results_directory/singlearea/output

echo submit-run-singlearea: compiling
Rscript epimap/compile.r

# echo submit-run-singlearea: running areas
# sbatch --wait \
#     --mail-user=$USER@stats.ox.ac.uk \
#     --mail-type=ALL \
#     --job-name=Rmap-singlearea \
#     --output=$results_directory/singlearea/output/run_%A_%a.out \
#     --partition=ziz-medium \
#     --ntasks=1 \
#     --time=18:00:00 \
#     --mem-per-cpu=5G \
#     --cpus-per-task=1 \
#     --array=1-$N \
#     --wrap \
#     "Rscript covidmap/stage1_run.r --area_index \$SLURM_ARRAY_TASK_ID $options"

echo submit-run-singlearea: combining areas
sbatch --wait \
    --mail-user=$USER@stats.ox.ac.uk \
    --mail-type=ALL \
    --job-name=Rmap-combineareas \
    --output=$results_directory/singlearea/output/combine_%A_%a.out \
    --partition=ziz-medium \
    --ntasks=1 \
    --time=18:00:00 \
    --mem-per-cpu=10G \
    --cpus-per-task=1 \
    --wrap \
    "Rscript covidmap/stage1_combine.r $options"

wait

echo submit-run-singlearea: DONE
