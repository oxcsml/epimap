# build the data
cd simulation/latent_epidemic
python pull_data.py
python tehtropolis_simulate.py
cd -

# clean data
./slurm/simulation/submit-clean.sh fits/simulation/clean simulation/latent_epidemic/tehtropolis 5

# run all models
./slurm/simulation/submit-run-cleaned-latents.sh            fits/simulation/cleaned-latent              fits/simulation/clean simulation/latent_epidemic/tehtropolis &
./slurm/simulation/submit-run-cleaned-latents-no-gp.sh      fits/simulation/cleaned-latent-no-gp        fits/simulation/clean simulation/latent_epidemic/tehtropolis &
./slurm/simulation/submit-run-cleaned-latents-no-metapop.sh fits/simulation/cleaned-latent-no-metapop   fits/simulation/clean simulation/latent_epidemic/tehtropolis &

./slurm/simulation/submit-run-latent-reports.sh             fits/simulation/latent-reports               fits/simulation/clean simulation/latent_epidemic/tehtropolis &
./slurm/simulation/submit-run-latent-reports-no-gp.sh       fits/simulation/latent-reports-no-gp         fits/simulation/clean simulation/latent_epidemic/tehtropolis &
./slurm/simulation/submit-run-latent-reports-no-metapop.sh  fits/simulation/latent-reports-no-metapop    fits/simulation/clean simulation/latent_epidemic/tehtropolis &

./slurm/simulation/submit-run-cori.sh                       fits/simulation/cori                        fits/simulation/clean simulation/latent_epidemic/tehtropolis &
