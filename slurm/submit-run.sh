#!/bin/bash

set -e

if [ $# -ne 1 ]; then
  echo Usage: submit-run results_directory
  exit 1
fi

results_directory=$1
echo results_directory = $results_directory

options="\
    --time_steps 15 \
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
    --mem-per-cpu=20G \
    --array=1-10 \
    --wrap \
    "Rscript mapping/run.r $options --cleaned_sample_id \$SLURM_ARRAY_TASK_ID"
    
wait

echo Combining results

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
    "Rscript mapping/merge_results.r $options"
wait

