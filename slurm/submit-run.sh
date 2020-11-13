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
    --fixed_gp_space_length_scale 3 \
    --fixed_gp_time_length_scale 21 \
"
slurm/submit-run-generic.sh $results_directory $options
