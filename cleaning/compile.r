library(rstan)

rstan_options(auto_write = TRUE)

stan_model(file = 'cleaning/stan_files/Rmap-clean.stan')
