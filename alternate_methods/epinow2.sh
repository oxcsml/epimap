#!/usr/bin/env bash

set -e

preprocess_pth=${code}/rmap/fits/epinow2/inputs.csv
outputs_folder=${code}/rmap/fits/epinow2/outputs
region_codes=${code}/rmap/fits/epinow2/region_codes.json
weeks_modelled=15
forecast_days=21
ncores=4 # number of cores for stan to use

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
        --output_folder=${folder} \
        --region_codes=${region_codes}

}

# This preprocessing step only needs to be done once, before all the runs
python epinow2_preprocess.py \
    ${code}/rmap/data/uk_cases.csv \
    ${preprocess_pth} \
    --region_codes=${region_codes}


# run the function single_run for each start date and area
start_date="2020-03-01" 
logsdir=${outputs_folder}/logs
mkdir -p $logsdir
logfile="${logsdir}/epinow2-logs-${start_date}.log"

for area in \
    "Adur" "Allerdale" "Amber Valley" "Arun" "Ashfield" "Ashford" "Babergh" "Barking and Dagenham" "Barnet" "Barnsley"
do
    echo "Running $area"

    echo "\n------------------------------\n" >> $logfile
    echo "\n------------------------------\n" >> $logfile
    echo "${area}" >> $logfile

    single_run ${start_date} "${area}" 2>&1 | tee -a $logfile > /dev/null
done

# post-processing to be applied to each date directory
python epinow2_postprocess.py \
     ${outputs_folder}/${start_date} \
     ${outputs_folder} \
     ${start_date} \
     ${forecast_days} \
     --region_codes=${region_codes} \
     --prefix=${start_date}

echo "Done"

