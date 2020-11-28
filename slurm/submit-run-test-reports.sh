#!/bin/bash
results_directory=fits/test-reports-1109
options=" \
    --fixed_gp_time_length_scale 100.0 \
    --fixed_gp_space_length_scale 3.0 \
    --clean_directory fits/clean-1109 \
    --observation_data latent_reports \
    --observation_model negative_binomial_2 \
"
    #--limit_area Oxford \
    #--limit_radius 0.8 \
slurm/submit-run-generic.sh $results_directory $options

