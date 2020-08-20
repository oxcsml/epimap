#!/usr/bin/bash

if [ "$#" -ne 2 ]; then
  echo reinflate.sh unique_experiment_id_in_fits experiment_id_in_website
else
  mkdir website/$2 
  python3 reinflate.py \
    `ls fits/*$1*_Rt.csv` `ls fits/*$1*_Cproj.csv` \
    website/$2/Rt.csv website/$2/Cproj.csv
  
fi

