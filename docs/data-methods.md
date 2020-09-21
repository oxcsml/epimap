---
layout: home
title: Data & Methods
weight: 1
---

We aim to develop a Bayesian method to reliably estimate Rt across local authorities with appropriately quantified uncertainties.

## Data

We use publicly available Pillar 1+2 case count time series:
*   312 lower-tier local authorities (LTLA) in England,
*   14 NHS Health Board level in Scotland, and
*   22 Unitary local authorities in Wales

## Model

At its core, our model uses a renewal equation formulation of epidemic dynamics within each local authority, building on the method of [Cori et al (2013)](https://doi.org/10.1093/aje/kwt133).


In addition,
*   We use a spatio-temporal Gaussian process to model correlations in effective reproduction numbers across neighbouring local authorities and across neighbouring points in time.
*   We use a cross-coupled metapopulation approach to model infections crossing local authority boundaries.

