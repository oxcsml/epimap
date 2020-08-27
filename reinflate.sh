#!/usr/bin/bash

if [ "$#" -ne 2 ]; then
  echo reinflate.sh unique_experiment_id_in_fits experiment_id_in_website
else
  mkdir website/$2 
  python3 reinflate.py \
    `ls fits/*$1*_Rt.csv` \
    `ls fits/*$1*_Pexceed.csv` \
    `ls fits/*$1*_Cpred.csv` \
    `ls fits/*$1*_Cproj.csv` \
    website/$2/Rt.csv \
    website/$2/Pexceed.csv \
    website/$2/Cpred.csv \
    website/$2/Cproj.csv

  cp -r website/$2 ~/pub_html/Rmap_time
  cp `ls fits/*$1*_pairs.pdf` ~/pub_html/Rmap_time/$2/pairs.pdf
  chmod a+rx ~/pub_html/Rmap_time/$2
  chmod a+r ~/pub_html/Rmap_time/$2/*
fi
