#!/bin/bash

trap 'echo submit-run-regional: Failed before finishing with exit code $? && exit $?' ERR

source ./slurm/cluster-config
echo "Compute cluster config: mail=$MAIL mail_type=$MAIL_TYPE partition=$PARTITION partition_large=$PARTITION_LARGE"

if [ $# == 1 ]
then
  results_directory=$1
  options="\
    --results_directory $results_directory \
  "
  N=9
elif [ $# == 2 ]
then
  results_directory=$1
  options="\
    --results_directory $results_directory \
    $2
  "
  N=9
elif [ $# == 3 ]
then
  results_directory=$1
  options="\
    --results_directory $results_directory \
    $2
  "
  N=$3
else
  echo Usage: submit-run-regional results_directory [options] [N]
  exit 1
fi

echo submit-run-regional: Inferring for each single area sample

echo results_directory = $results_directory
mkdir -p $results_directory
mkdir -p $results_directory/regional
mkdir -p $results_directory/regional/output
git rev-parse HEAD > $results_directory/regional/git-hash.txt
options="--approximation regional --num_regions $N $options"
echo $options

echo submit-run-regional: compiling
Rscript epimap/compile.r

echo submit-run-regional: running regions
sbatch --wait \
    --mail-user=$MAIL \
    --mail-type=$MAIL_TYPE \
    --job-name=Rmap-regional \
    --output=$results_directory/regional/output/run_%A_%a.out \
    --partition=$PARTITION_LARGE \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=10G \
    --array=1-$N \
    --wrap \
    "Rscript covidmap/stage2_run.r --region_id \$SLURM_ARRAY_TASK_ID $options"

echo submit-run-regional: Merging results

sbatch --wait \
    --mail-user=$MAIL \
    --mail-type=$MAIL_TYPE \
    --job-name=Rmap-mergeregions \
    --output=$results_directory/regional/output/merge_%A_%a.out \
    --partition=$PARTITION_LARGE \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=20G \
    --wrap \
    "Rscript covidmap/stage2_merge.r $options"

echo submit-run-regional: DONE
