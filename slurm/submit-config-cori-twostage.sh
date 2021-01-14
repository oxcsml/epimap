#!/bin/bash

trap 'echo submit-run: Failed before finishing with exit code $? && exit $?' ERR

if [ $# -ne 2 ]; then
  echo Usage: submit-run results_directory singlearea_directory
  exit 1
fi

results_directory=$1
singlearea_directory=$2

options="\
    --time_steps 15 \
    --iterations 6000 \
    --observation_data cleaned_recon_sample \
    --observation_model poisson \
    --spatialkernel none \
    --temporalkernel none \
    --localkernel local \
    --globalkernel none \
    --metapop none"

./slurm/submit-run-twostage.sh $results_directory $singlearea_directory $options