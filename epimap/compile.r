library(rstan)

rstan_options(auto_write = TRUE)

stan_model(file = 'epimap/stan_files/Rmap-clean.stan')
stan_model(file = 'epimap/stan_files/Rmap.stan')
stan_model(file = 'epimap/stan_files/Rmap-latent.stan')
