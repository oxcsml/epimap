---
layout: home
title: Data & Methods
weight: 1
---

We aim to develop a Bayesian method to reliably estimate $R_t$ across local authorities in the UK with appropriately quantified uncertainties.

## Data

We use publicly available Pillar 1+2 case count time series:
*   312 lower-tier local authorities (LTLA) in England,
*   14 NHS Health Board level in Scotland, and
*   22 Unitary local authorities in Wales

## Methods

At its core, our model uses a renewal equation formulation of epidemic dynamics within each local authority, building on the method of [Cori et al (2013)](https://doi.org/10.1093/aje/kwt133).


In addition,
*   Correlations in effective reproduction numbers across neighbouring local authorities and across neighbouring points in time are modelled using a spatio-temporal Gaussian process.
*   Potential infections that cross local authority boundaries are accounted for using a cross-coupled metapopulation approach. In order to do so, we incorporate both real commuter data as well as a well-known [human mobility model](https://arxiv.org/abs/1111.0586).
*   Problems associated uncertainties in case report, outliers in case counts and delays in the testing procedure are alleviated via an autoregressive latent process that preprocesses observed case counts.
*   We use the [No-U-Turn Sampler](https://arxiv.org/abs/1111.4246) within the [Stan](https://mc-stan.org/) probabilistic programming language for posterior inference for both our model and preprocessing

