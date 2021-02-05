#!/bin/bash

trap 'echo submit-run-regional: Failed before finishing with exit code $? && exit $?' ERR

if [ $# -lt 2 ]; then
  echo Usage: submit-run-regional results_directory options
  exit 1
fi

echo submit-run-regional: Inferring for each single area sample

results_directory=$1
echo results_directory = $results_directory
mkdir -p $results_directory
mkdir -p $results_directory/regional
mkdir -p $results_directory/regional/output
git rev-parse HEAD > $results_directory/regional/git-hash.txt
options="--approximation regional --results_directory $results_directory ${@:2}"
echo $options

echo submit-run-regional: compiling
Rscript epimap/compile.r

echo submit-run-regional: running regions
sbatch --wait \
    --mail-user=$USER@stats.ox.ac.uk \
    --mail-type=ALL \
    --job-name=Rmap-regional \
    --output=$results_directory/regional/output/run_%A_%a.out \
    --partition=ziz-large \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=10G \
    --array=1-9 \
    --wrap \
    "Rscript covidmap/stage2_run.r $options --region_id \$SLURM_ARRAY_TASK_ID && echo Rmap-regional: DONE"

echo submit-run-regional: Merging results

sbatch --wait \
    --mail-user=$USER@stats.ox.ac.uk \
    --mail-type=ALL \
    --job-name=Rmap-mergeregions \
    --output=$results_directory/regional/output/merge_%A_%a.out \
    --partition=ziz-large \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=20G \
    --wrap \
    "Rscript covidmap/stage2_merge.r $options && echo Rmap-mergeregions: DONE"

echo submit-run-regional: DONE
