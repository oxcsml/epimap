#!/bin/bash

jobname=$(date +'%Y-%m-%d')
results_directory="fits/Rmap-$jobname"
echo results_directory = $results_directory

options="\
    --time_steps 25 \
    --iterations 6000 \
    --observation cleaned_recon_sample \
    --results_directory $results_directory"

sbatch --wait \
    --mail-user=$USER@stats.ox.ac.uk \
    --mail-type=ALL \
    --job-name=Rmap_run \
    --output=slurm/output/run_%A_%a.out \
    --partition=ziz-large \
    --ntasks=1 \
    --time=40:00:00 \
    --mem-per-cpu=30G \
    --cpus-per-task=1 \
    --array=1-10 \
    --wrap \
    "Rscript run.r $options --cleaned_sample_id \$SLURM_ARRAY_TASK_ID"
wait

echo Combining results

Rscript postprocess_samples.r $options

mkdir -p docs/assets/data/${jobname}

python3 reinflate.py \
    ${results_directory}/merged_Rt.csv \
    ${results_directory}/merged_Pexceed.csv \
    ${results_directory}/merged_Cweekly.csv \
    ${results_directory}/merged_Cpred.csv \
    ${results_directory}/merged_Cproj.csv \
    docs/assets/data/${jobname}/Rt.csv \
    docs/assets/data/${jobname}/Pexceed.csv \
    docs/assets/data/${jobname}/Cweekly.csv \
    docs/assets/data/${jobname}/Cpred.csv \
    docs/assets/data/${jobname}/Cproj.csv

echo Done
