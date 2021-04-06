#!/bin/bash

trap 'echo submit-twostage: Failed before finishing with exit code $? && exit $?' ERR

if [ $# == 1 ]
then
  results_directory=$1
  options="\
    --results_directory $results_directory \
  "
  N=9
elif [ $# == 2 ]
then
  results_directory=$1
  options="\
    --results_directory $results_directory \
    $2
  "
  N=9
elif [ $# == 3 ]
then
  results_directory=$1
  options="\
    --results_directory $results_directory \
    $2
  "
  N=$3
else
  echo Usage: submit-run-regional results_directory [options] [N]
  exit 1
fi

echo submit-run: Inferring for each single area sample

results_directory=$1
echo results_directory = $results_directory
mkdir -p $results_directory
mkdir -p $results_directory/twostage
mkdir -p $results_directory/twostage/output
git rev-parse HEAD > $results_directory/git-hash.txt
options="--approximation twostage --num_samples $N $options"
echo $options

echo submit-run-twostage: compiling
Rscript epimap/compile.r

echo submit-run-twostage: running samples
sbatch --wait \
    --mail-user=$USER@stats.ox.ac.uk \
    --mail-type=ALL \
    --job-name=Rmap_run \
    --output=$results_directory/twostage/output/run_%A_%a.out \
    --partition=ziz-large \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=10G \
    --array=1-$N \
    --wrap \
    "Rscript covidmap/stage2_run.r \
        --singlearea_sample_id \$SLURM_ARRAY_TASK_ID \
        $options"

echo submit-run-twostage: Merging results

sbatch --wait \
    --mail-user=$USER@stats.ox.ac.uk \
    --mail-type=ALL \
    --job-name=Rmap_merge \
    --output=$results_directory/twostage/output/merge_%A_%a.out \
    --partition=ziz-large \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=20G \
    --wrap \
    "Rscript covidmap/stage2_merge.r $options"

echo submit-run: DONE
