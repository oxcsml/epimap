#!/bin/bash

trap 'echo submit-run: Failed before finishing with exit code $? && exit $?' ERR

if [ $# -ne 2 ]; then
  echo Usage: submit-run results_directory singlearea_directory
  exit 1
fi

results_directory=$1
singlearea_directory=$2

options="\
    --iterations 6000 \
    --observation_data latent_reports \
    --observation_model gaussian \
"
./slurm/submit-run-twostage.sh $results_directory $singlearea_directory $options