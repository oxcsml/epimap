#!/bin/bash

trap 'echo submit-run: Failed before finishing with exit code $? && exit $?' ERR

if [ $# -ne 3 ]; then
  echo Usage: submit-run results_directory clean_directory data_directory
  exit 1
fi

results_directory=$1
clean_directory=$2
data_directory=$3
options="\
    --weeks_modelled 29 \
    --iterations 6000 \
    --observation_data cleaned_recon_sample \
    --observation_model negative_binomial_3 \
    --clean_directory $clean_directory \
    --data_directory $data_directory \
"
slurm/submit-run-generic.sh $results_directory $options
