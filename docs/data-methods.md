---
layout: page
title: Data & Methods
weight: 1
---

We aim to reliably estimate the instantaneous reproduction number $R_t$ across local authorities in the UK with appropriately quantified uncertainties.

### Data

We use publicly available [Pillar 1+2](https://www.gov.uk/government/publications/coronavirus-covid-19-scaling-up-testing-programmes) 
daily counts of positive PCR swab tests by specimen date, for:
*   312 lower-tier local authorities (LTLA) in England 
([here](https://coronavirus.data.gov.uk)),
*   14 NHS Health Board level in Scotland (each covering multiple local authorities)
([here](https://www.gov.scot/publications/coronavirus-covid-19-daily-data-for-scotland/)), and
*   22 Unitary local authorities in Wales ([here](https://phw.nhs.wales/topics/latest-information-on-novel-coronavirus-covid-19/)).

Other data sources:
* UK 2011 Census commuter flow data ([here](https://www.statistics.digitalresources.jisc.ac.uk)),
* ONS UK population estimates from mid 2019 ([here](https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesforukenglandandwalesscotlandandnorthernireland)).


### Methods

At its core, our Bayesian method uses a renewal equation formulation of epidemic dynamics 
within each local authority, 
building on the methods of [Cori et al (2013)](https://doi.org/10.1093/aje/kwt133) and 
[Flaxman et al (2020)](https://www.nature.com/articles/s41586-020-2405-7).

Specific extensions that we have employed to adapt the renewal equation approach to our local-level model:
*   Correlations in effective $R_t$ across neighbouring local authorities and across neighbouring points in time are modelled using a spatio-temporal Gaussian process. This allows for sharing of statistical strengths.
*   Potential infections that cross local authority boundaries are accounted for using a cross-coupled metapopulation approach. In order to do so, we incorporate real commuter data from the UK 2011 Census ([here](https://www.statistics.digitalresources.jisc.ac.uk)).
*   Problems associated with noise in the case reporting process, outliers in case counts and delays in the testing and reporting system are alleviated by modelling the epidemic using a latent process with associated observation model for reported cases, following [Flaxman et al (2020)](https://www.nature.com/articles/s41586-020-2405-7).
*   We use the [No-U-Turn Sampler](https://arxiv.org/abs/1111.4246) inplemented in the [Stan](https://mc-stan.org/) probabilistic programming system for posterior inference.
*   Because of the computational cost of posterior simulation in the resulting complex model, 
we split the posterior simulation into two phases: 
an initial phase which infers the latent epidemic process in each local authority, 
and a second phase which infers the $R_t$ and metapopulation model parameters.

Detailed description of the method and source code will be provided as soon as possible.


