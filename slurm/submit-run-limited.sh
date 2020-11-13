#!/bin/bash

trap 'echo submit-run: Failed before finishing with exit code $? && exit $?' ERR

if [ $# -ne 2 ]; then
  echo Usage: submit-run results_directory clean_directory
  exit 1
fi

results_directory=$1
clean_directory=$2
options="\
    --weeks_modelled 15 \
    --iterations 6000 \
    --observation_data cleaned_latent_sample \
    --observation_model gaussian \
    --clean_directory $clean_directory \
    --limit_area Oxford \
    --limit_radius 0.8 \
    --fixed_gp_space_length_scale 3 \
    --fixed_gp_time_length_scale 21 \
    --days_predicted 3 \
    --num_steps_forecasted 4 \
"
slurm/submit-run-generic.sh $results_directory $options


