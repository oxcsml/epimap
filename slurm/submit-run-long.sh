#!/bin/bash

trap 'echo submit-run: Failed before finishing with exit code $? && exit $?' ERR

if [ $# -ne 1 ]; then
  echo Usage: submit-run-long results_directory
  exit 1
fi

echo submit-run: Inferring for each cleaned sample

results_directory=$1
echo results_directory = $results_directory

options="\
    --time_steps 28 \
    --iterations 8000 \
    --observation cleaned_recon_sample \
    --results_directory $results_directory"

sbatch --wait \
    --mail-user=$USER@stats.ox.ac.uk \
    --mail-type=ALL \
    --job-name=Rmap_run \
    --output=slurm/output/run_%A_%a.out \
    --partition=ziz-large \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=40G \
    --array=1-10 \
    --wrap \
    "Rscript mapping/run.r $options \
        --cleaned_sample_id \$SLURM_ARRAY_TASK_ID \
        && echo run: DONE"

echo submit-run: Merging results

sbatch --wait \
    --mail-user=$USER@stats.ox.ac.uk \
    --mail-type=ALL \
    --job-name=Rmap_merge \
    --output=slurm/output/merge_%A_%a.out \
    --partition=ziz-large \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=30G \
    --wrap \
    "Rscript mapping/merge_results.r $options && echo merge: DONE"

echo submit-run: DONE
