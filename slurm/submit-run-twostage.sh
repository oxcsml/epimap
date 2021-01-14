#!/bin/bash

trap 'echo submit-twostage: Failed before finishing with exit code $? && exit $?' ERR

if [ $# -lt 2 ]; then
  echo Usage: submit-run-twostage results_directory singlearea_directory options
  exit 1
fi

echo submit-run: Inferring for each single area sample

results_directory=$1
echo results_directory = $results_directory
mkdir -p $results_directory
mkdir -p $results_directory/output
git rev-parse HEAD > $results_directory/git-hash.txt
options="--approximation twostage --results_directory $results_directory --singlearea_directory $2 ${@:3}"
echo $options

echo submit-run-twostage: compiling
Rscript mapping/compile.r

echo submit-run-twostage: running samples
sbatch --wait \
    --mail-user=$USER@stats.ox.ac.uk \
    --mail-type=ALL \
    --job-name=Rmap_run \
    --output=$results_directory/output/run_%A_%a.out \
    --partition=ziz-large \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=10G \
    --array=1-10 \
    --wrap \
    "Rscript covidmap/stage2_run.r $options \
        --singlearea_sample_id \$SLURM_ARRAY_TASK_ID \
        && echo run: DONE"

echo submit-run-twostage: Merging results

sbatch --wait \
    --mail-user=$USER@stats.ox.ac.uk \
    --mail-type=ALL \
    --job-name=Rmap_merge \
    --output=$results_directory/output/merge_%A_%a.out \
    --partition=ziz-large \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=20G \
    --wrap \
    "Rscript covidmap/stage2_merge.r $options && echo merge: DONE"

echo submit-run: DONE
