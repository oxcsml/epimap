#!/bin/bash

numiters=4000
numchains=4

Rscript run.r -i $numiters -c $numchains -k none -m none -l poisson
Rscript run.r -i $numiters -c $numchains -k none -m none -l negative_binomial
Rscript run.r -i $numiters -c $numchains -k matern12 -m none -l negative_binomial
Rscript run.r -i $numiters -c $numchains -k none -m uniform1 -l negative_binomial
Rscript run.r -i $numiters -c $numchains -k none -m uniform2 -l negative_binomial
Rscript run.r -i $numiters -c $numchains -k matern32 -m none -l negative_binomial
Rscript run.r -i $numiters -c $numchains -k exp_quad -m none -l negative_binomial
Rscript run.r -i $numiters -c $numchains -k matern12 -m uniform1 -l negative_binomial
Rscript run.r -i $numiters -c $numchains -k matern12 -m uniform2 -l negative_binomial
