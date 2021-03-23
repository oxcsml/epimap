mkdir -p simulation_fits/test

cp simulation/latent_epidemic/tehtropolis/sample/cases.csv simulation_fits/test

./slurm/submit-run-epiestim.sh \
    simulation_fits/test \
    "--data_directory simulation/latent_epidemic/tehtropolis/sample \
    --first_day_modelled 2020-04-09 \
    --days_ignored 21" 5

./slurm/submit-run-singlearea.sh \
    simulation_fits/test \
    "--data_directory simulation/latent_epidemic/tehtropolis/sample \
    --first_day_modelled 2020-04-09 \
    --days_ignored 21 \
    --days_predicted 21 \
    --num_steps_forecasted 3" 5

./slurm/submit-run-twostage.sh \
    simulation_fits/test \
    "--data_directory simulation/latent_epidemic/tehtropolis/sample \
    --first_day_modelled 2020-04-09 \
    --days_ignored 21 \
    --days_predicted 21 \
    --num_steps_forecasted 3" 10

./slurm/submit-run-regional.sh \
    simulation_fits/test \
    "--data_directory simulation/latent_epidemic/tehtropolis/sample \
    --first_day_modelled 2020-04-09 \
    --days_ignored 21 \
    --days_predicted 21 \
    --num_steps_forecasted 3" 1