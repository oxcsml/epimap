#!/bin/bash
trap 'echo submit-bootstrap-runs-regional: Failed before finishing with exit code $? && exit $?' ERR

# Activate the right bash environment

echo $#
if [ $# == 3 ]
then
    N_bootstrap=$1
    results_directory=$2
    options=$3
else
  echo Usage: submit-runs-bootstrap.sh num_bootstrap results_directory options
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
    echo "Running stage 2 regional for bootstrap sample ${i}"
    slurm/submit-run-regional.sh $sample_results_directory "${sample_options}" &
    sleep 15
done
wait

# MERGE
mkdir -p $results_directory/regional/output

N_regions=9
options="\
    --results_directory $results_directory \
    --approximation regional --num_regions $N_regions \
    --num_samples $N_bootstrap \
    $options
"
options="--bootstrap_merge TRUE $options"

echo daily_bootstrap_update: merging bootstrap samples
sbatch --wait \
    --mail-user=$USER@stats.ox.ac.uk \
    --mail-type=ALL \
    --job-name=Rmap-mergeregions_bootstrap \
    --output=$results_directory/regional/output/merge_%A_%a.out \
    --partition=ziz-large \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=20G \
    --wrap \
    "Rscript covidmap/stage2_merge.r ${options}"

echo completed
