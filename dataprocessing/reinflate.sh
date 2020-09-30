#!/bin/bash

set -e

if [ $# -ne 2 ]; then
  echo Usage: reinflate results-directory website-folder
  exit 1
fi

results_directory=$1
website_folder=$2

mkdir -p docs/assets/data/${website_folder}

python3 data/scripts/reinflate.py \
    ${results_directory}/merged_Rt.csv \
    ${results_directory}/merged_Pexceed.csv \
    ${results_directory}/merged_Cweekly.csv \
    ${results_directory}/merged_Cpred.csv \
    ${results_directory}/merged_Cproj.csv \
    docs/assets/data/${website_folder}/Rt.csv \
    docs/assets/data/${website_folder}/Pexceed.csv \
    docs/assets/data/${website_folder}/Cweekly.csv \
    docs/assets/data/${website_folder}/Cpred.csv \
    docs/assets/data/${website_folder}/Cproj.csv

echo Done