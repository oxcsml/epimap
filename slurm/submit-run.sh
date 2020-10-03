#!/bin/bash

trap 'echo submit-run: Failed before finishing with exit code $? && exit $?' ERR

if [ $# -ne 2 ]; then
  echo Usage: submit-run results_directory clean_directory
  exit 1
fi

echo submit-run: Inferring for each cleaned sample

results_directory=$1
clean_directory=$2
echo results_directory = $results_directory
echo clean_directory = $clean_directory

mkdir -p $results_directory/output

options="\
    --time_steps 15 \
    --iterations 6000 \
    --observation_data cleaned_recon_sample \
    --observation_model negative_binomial_3 \
    --results_directory $results_directory \
    --clean_directory $clean_directory"

sbatch --wait \
    --mail-user=$USER@stats.ox.ac.uk \
    --mail-type=ALL \
    --job-name=Rmap_run \
    --output=$results_directory/output/run_%A_%a.out \
    --partition=ziz-large \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=20G \
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
    --output=$results_directory/output/merge_%A_%a.out \
    --partition=ziz-large \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=100G \
    --wrap \
    "Rscript mapping/merge_results.r $options && echo merge: DONE"

echo submit-run: DONE
