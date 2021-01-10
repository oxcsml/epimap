#!/bin/bash

trap 'echo submit-run-full: Failed before finishing with exit code $? && exit $?' ERR

if [ $# -lt 2 ]; then
  echo Usage: submit-run-full results_directory options
  exit 1
fi

echo submit-run: Inferring for each cleaned sample

results_directory=$1
echo results_directory = $results_directory
mkdir -p $results_directory
mkdir -p $results_directory/output
git rev-parse HEAD > $results_directory/git-hash.txt
options="--stage full --results_directory $results_directory ${@:2}"
echo $options

echo compiling
Rscript mapping/compile.r

sbatch --wait \
    --mail-user=$USER@stats.ox.ac.uk \
    --mail-type=ALL \
    --job-name=Rmap_run \
    --output=$results_directory/output/run_%A_%a.out \
    --partition=ziz-large \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=10G \
    --array=1-9 \
    --wrap \
    "Rscript mapping/run.r $options \
        --region_id \$SLURM_ARRAY_TASK_ID \
        && echo run: DONE"

echo submit-run: Merging results

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
    "Rscript mapping/merge_results.r $options && echo merge: DONE"

echo submit-run: DONE
