#!/usr/bin/bash

if [ "$#" -ne 2 ]; then
  echo reinflate.sh unique_experiment_id_in_fits experiment_id_in_docs/assets/data/
else
  mkdir docs/assets/data/$2 
  python3 reinflate.py \
    `ls fits/*$1*_Rt.csv` \
    `ls fits/*$1*_Pexceed.csv` \
    `ls fits/*$1*_Cweekly.csv` \
    `ls fits/*$1*_Cpred.csv` \
    `ls fits/*$1*_Cproj.csv` \
    docs/assets/data/$2/Rt.csv \
    docs/assets/data/$2/Pexceed.csv \
    docs/assets/data/$2/Cweekly.csv \
    docs/assets/data/$2/Cpred.csv \
    docs/assets/data/$2/Cproj.csv

  cp -r docs/assets/data/$2 ~/pub_html/Rmap
  #cp `ls fits/*$1*_pairs.pdf` ~/pub_html/Rmap/$2/pairs.pdf
  chmod a+rx ~/pub_html/Rmap/$2
  chmod a+r ~/pub_html/Rmap/$2/*
fi
