library(rstan)

rstan_options(auto_write = TRUE)

stan_model(file = 'mapping/stan_files/Rmap.stan')
stan_model(file = 'mapping/stan_files/Rmap-latent.stan')
