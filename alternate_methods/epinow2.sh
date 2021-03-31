
preprocess_pth=${code}/rmap/fits/epinow2/inputs.csv
outputs_folder=${code}/rmap/fits/epinow2/outputs
weeks_modelled=15
forecast_days=21
ncores=4 # number of cores for stan to use

# This preprocessing step only needs to be done once, before all the runs
python epinow2_preprocess.py ${code}/rmap/data/uk_cases.csv ${preprocess_pth} 

echo "preprocessing done"

single_run() { # $1=start_date, $2=area

    folder="${outputs_folder}/$1"
    mkdir -p $folder

    Rscript epinow2_run.r \
        --first_day_modelled="$1" \
        --weeks_modelled=${weeks_modelled} \
        --area="$2" \
        --case_counts=${preprocess_pth} \
        --ncores=${ncores} \
        --forecast_horizon=${forecast_days} \
        --output_folder=${folder}
}

# do this for each start date, area
single_run "2021-03-01" "Adur"

# post-processing to be applied to each directory

echo "Done"

