#!/bin/bash
results_directory=fits/test-reports-fixed
options=" \
    --fixed_gp_time_length_scale 60.0 \
    --fixed_gp_space_length_scale 3.0 \
    --clean_directory fits/clean-weekly \
    --observation_data latent_reports \
    --observation_model negative_binomial_2 \
"
slurm/submit-run-generic.sh $results_directory $options

