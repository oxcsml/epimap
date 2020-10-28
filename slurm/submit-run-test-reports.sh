#!/bin/bash
results_directory=fits/test-reports
options=" \
    --clean_directory fits/clean-newprior \
    --observation_data latent_reports \
    --observation_model negative_binomial_2 \
"
slurm/submit-run-generic.sh $results_directory $options

