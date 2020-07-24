#!/bin/bash

numiters=4000
numchains=4

# Rscript run.r -i $numiters -c $numchains -s none -m none -o poisson
# Rscript run.r -i $numiters -c $numchains -s none -m none -o negative_binomial
# Rscript run.r -i $numiters -c $numchains -s matern12 -m none -o negative_binomial
# Rscript run.r -i $numiters -c $numchains -s none -m uniform1 -o negative_binomial
# Rscript run.r -i $numiters -c $numchains -s none -m uniform2 -o negative_binomial
# Rscript run.r -i $numiters -c $numchains -s matern32 -m none -o negative_binomial
# Rscript run.r -i $numiters -c $numchains -s exp_quad -m none -o negative_binomial
# Rscript run.r -i $numiters -c $numchains -s matern12 -m uniform1 -o negative_binomial
# Rscript run.r -i $numiters -c $numchains -s matern12 -m uniform2 -o negative_binomial

Rscript run.r -i $numiters -c $numchains -s matern12 -l local -m uniform1 -o negative_binomial
Rscript run.r -i $numiters -c $numchains -s matern12 -l none -m uniform1 -o negative_binomial

