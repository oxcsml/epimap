#!/bin/bash

trap 'echo submit-run: Failed before finishing with exit code $? && exit $?' ERR

if [ $# -le 2 ]; then
  echo Usage: submit-run results_directory clean_directory
  exit 1
fi

results_directory=$1
clean_directory=$2
options="\
    --clean_directory $clean_directory \
"
slurm/submit-run-generic.sh $results_directory $options
