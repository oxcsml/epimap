#!/bin/bash
results_directory=fits/test-latents-fixed-2
options=" \
    --fixed_gp_time_length_scale 60.0 \
    --fixed_gp_space_length_scale 3.0 \
    --clean_directory fits/clean-weekly \
    --observation_data cleaned_latent_sample \
    --observation_model gaussian \
"
slurm/submit-run-generic.sh $results_directory $options

