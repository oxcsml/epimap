sbatch --wait \
    --mail-user=$USER@stats.ox.ac.uk \
    --mail-type=ALL \
    --job-name=clean_ts \
    --output=slurm/output/cleaning/clean_timeseries_%A_%a.out \
    --partition=ziz-medium \
    --ntasks=1 \
    --time=18:00:00 \
    --mem-per-cpu=5G \
    --cpus-per-task=1 \
    --array=1-348 \
    --wrap \
    'Rscript clean_area.r --task_id $SLURM_ARRAY_TASK_ID'
wait

echo Combining results...

Rscript clean_combine.r
