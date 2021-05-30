#!/bin/bash
trap 'echo submit-bootstrap-runs-singlearea: Failed before finishing with exit code $? && exit $?' ERR

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
    echo "Running stage 1 singlearea for bootstrap sample ${i}"
    slurm/submit-run-singlearea.sh $sample_results_directory "${sample_options}" &
    sleep 15
done
wait

# MERGE
options="\
    --produce_plots FALSE \
    --results_directory $results_directory \
    $options
"
options="--num_bootstrap $N_bootstrap $options"

mkdir -p $results_directory/singlearea/output

echo daily_bootstrap_update: merging bootstrap samples
sbatch --mail-user=$USER@stats.ox.ac.uk \
    --mail-type=ALL \
    --job-name=Rmap-merge-singlearea \
    --output=$results_directory/singlearea/output/merge_%A_%a.out \
    --partition=ziz-large \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=20G \
    --wrap \
    "Rscript covidmap/stage1_bootstrap_combine.r ${options}"
echo completed
