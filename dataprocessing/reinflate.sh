#!/bin/bash

set -e

if [ $# -ne 2 ]; then
  echo Usage: reinflate results-prefix website-folder
  exit 1
fi

results_prefix=$1
website_folder=$2

mkdir -p docs/assets/data/${website_folder}

python3 dataprocessing/reinflate.py \
    ${results_prefix}Rt.csv \
    ${results_prefix}Pexceed.csv \
    ${results_prefix}Bweekly.csv \
    ${results_prefix}Cweekly.csv \
    ${results_prefix}Bpred.csv \
    ${results_prefix}Bproj.csv \
    ${results_prefix}Cpred.csv \
    ${results_prefix}Cproj.csv \
    ${results_prefix}Xpred.csv \
    ${results_prefix}Xproj.csv \
    docs/assets/data/${website_folder}/Rt.csv \
    docs/assets/data/${website_folder}/Pexceed.csv \
    docs/assets/data/${website_folder}/Bweekly.csv \
    docs/assets/data/${website_folder}/Cweekly.csv \
    docs/assets/data/${website_folder}/Bpred.csv \
    docs/assets/data/${website_folder}/Bproj.csv \
    docs/assets/data/${website_folder}/Cpred.csv \
    docs/assets/data/${website_folder}/Cproj.csv \
    docs/assets/data/${website_folder}/Xpred.csv \
    docs/assets/data/${website_folder}/Xproj.csv

echo reinflate: DONE
