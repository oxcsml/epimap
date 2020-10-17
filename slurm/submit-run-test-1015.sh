#!/bin/bash
results_directory=fits/test-1015
options=" \
    --days_ignored 0 \
    --days_predicted 1 \
    --iterations 3000 \
    --clean_directory fits/clean-reflectnormal \
    --observation_data cleaned_latent_sample \
    --observation_model gaussian \
"
slurm/submit-run-generic.sh $results_directory $options

