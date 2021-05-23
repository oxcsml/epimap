#!/bin/bash
trap 'echo submit-runs-bootstrap: Failed before finishing with exit code $? && exit $?' ERR

# Activate the right bash environment

echo $#
if [ $# == 4 ]
then
    N_bootstrap=$1
    results_directory=$2
    options=$3
    stage=$4
else
  echo Usage: submit-runs-bootstrap.sh num_bootstrap results_directory options stage
  exit 1
fi
for ((i=1; i <= ${N_bootstrap}; i++));
do  
    sample_results_directory="${results_directory}/bootstrap_$i"  
    mkdir -p $sample_results_directory
    cp data/cases.csv $sample_results_directory
    sample_options="\
        $(sed -n "${i}p" < ${results_directory}/bootstrap_params.txt) \
        ${options}
    "
    echo $sample_options
    if [ $stage == 1 ]; then
        echo "Running stage 1 singlearea for bootstrap sample ${i}"
        slurm/submit-run-singlearea.sh $sample_results_directory "${sample_options}" &
    elif [ $stage == 2 ]; then
        echo "Running stage 2 regional for bootstrap sample ${i}"
        slurm/submit-run-regional.sh $sample_results_directory "${sample_options}" &
    else
        echo "Stage ${stage} not found."
    fi
    sleep 15
done
wait
echo completed
