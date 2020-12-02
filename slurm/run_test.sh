Rscript mapping/compile.r

clean_directory=fits/clean-2020-10-25
results_directory=fits/testing/output

options="\
    --time_steps 15 \
    --iterations 6000 \
    --observation_data cleaned_latent_sample \
    --observation_model gaussian \
    --clean_directory $clean_directory \
    --results_directory $results_directory \
    --cleaned_sample_id 1 \
    --fixed_gp_time_length_scale 28 \
    --days_predicted 0 \
    --num_steps_forecasted 1 \
"

Rscript mapping/run.r $options

echo Done