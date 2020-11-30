#!/bin/bash

trap 'echo submit-run: Failed before finishing with exit code $? && exit $?' ERR

if [ $# -ne 2 ]; then
  echo Usage: submit-run results_directory clean_directory
  exit 1
fi

results_directory=$1
clean_directory=$2
options="\
    --iterations 6000 \
    --observation_data latent_reports \
    --observation_model gaussian \
    --clean_directory $clean_directory \
    --limit_area Oxford \
    --limit_radius 0.5 \
"
slurm/submit-run-generic.sh $results_directory $options


