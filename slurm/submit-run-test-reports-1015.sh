#!/bin/bash
results_directory=fits/test-reports-1015
options=" \
    --days_ignored 0 \
    --days_predicted 1 \
    --iterations 3000 \
    --clean_directory fits/clean-reflectnormal \
    --observation_data latent_reports \
    --observation_model negative_binomial_3 \
"
slurm/submit-run-generic.sh $results_directory $options

