#!/bin/bash
results_directory=fits/test-latents-1109
options=" \
    --first_day_modelled 2020-08-01 \
    --last_day_modelled 2020-10-30 \
    --iterations 3000 \
    --num_steps_forecasted 3 \
    --fixed_gp_time_length_scale 100.0 \
    --fixed_gp_space_length_scale 3.0 \
    --clean_directory fits/clean-1109 \
    --observation_data cleaned_latent_sample \
    --observation_model gaussian \
"
    #--limit_area Oxford \
    #--limit_radius 0.8 \
slurm/submit-run-generic.sh $results_directory $options

