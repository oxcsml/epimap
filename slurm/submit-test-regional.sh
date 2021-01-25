#!/bin/bash
results_directory=fits/test
options=" \
    --iterations 30 \
    --globalkernel none \
    --spatialkernel matern12 \
    --fixed_gp_time_length_scale 100.0 \
    --fixed_gp_space_length_scale 0.1 \
"
    #--first_day_modelled 2020-08-01 \
    #--last_day_modelled 2020-10-30 \
    #--limit_area Oxford \
    #--limit_radius 0.8 \
slurm/submit-run-regional.sh $results_directory $options

