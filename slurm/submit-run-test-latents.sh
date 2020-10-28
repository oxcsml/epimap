#!/bin/bash
results_directory=fits/test-latents
options=" \
    --clean_directory fits/clean-newprior \
    --observation_data cleaned_latent_sample \
    --observation_model gaussian \
"
slurm/submit-run-generic.sh $results_directory $options

