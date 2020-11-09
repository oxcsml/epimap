library(rstan)
library(gsubfn)
library(plyr)

#' Epiclean inference of Rt and incidence curve.
#'
#' Uses Bayesian inference, with a latent incidence process given by a 
#' continuous approximation of a negative-binomial renewal process, 
#' a piecewise constant AR1 prior on Rt, and an observation model for
#' counts of positive diagnoses described by a negative-binomial and a
#' delay distribution.
#' 
#' Implemented as MCMC sampling (No-U-Turns Sampler) in stan.
#'
#' @param Count Dataframe or list of case counts of positive tests by day; colnames(Count) can be date strings in YYYY-MM-DD format.
#' @param Tcond Initial number of days that are not modelled (but are conditioned on)
#' @param Nstep Number of AR1 time steps to model 
#' @param Nproj Number of AR1 time steps to forecast
#' @param Tstep Number of days per AR1 time step; default 1
#' @param infprofile Infection profile (generation interval distribution)
#' @param testdelayprofile Distribution of delay between infection and getting tested
#' @param resultdelayalpha Parameters of Dirichlet prior of delay distribution between getting tested and resulted being reported
#' @param gp_time_scale Length scale of AR1 process; default 14
#' @param gp_time_decay_scale Decay scale of AR1 process; default .1
#' @param fixed_gp_time_length_scale If specified, fixes AR1 length scale as specified; not specified by default
#' @param mu_scale Scale of prior on the prior mean of log(Rt)
#' @param sigma_scale Scale of prior on the prior standard deviation of log(Rt)
#' @param phi_latent_scale Scale of prior on negative-binomial dispersion parameter of latent epidemic process; default 5.0
#' @param phi_observed_scale Scale of prior on dispersion parameter of observation process of positive test counts; default 5.0
#' @param xi_scale Scale of prior on exegeneous infection rate
#' @param reconstruct_infections Whether to reconstruct infection day for each case; default True
#' @param outlier_prob_threshold Probability hreshold to determine if a diagnosis count is an outlier; default 1.0
#' @param outlier_count_threshold Count threshold to determine if a diagnosis count is an outlier; default 10
#' @param num_iterations Number of MCMC iterations; default 3000
#' @param num_chains Number of MCMC chains; default 1
#' @param percentiles Thresholds of percentiles for posterior statistics; default c(.025,.25,.5,.75,.975)
#'
#' @return A list consists of the stanfit MCMC output, as well as posterior 
#' statistics of the instantaneous reproduction numbers Rt and 
#' incidence numbers Xt.
#'
#' @export
epiclean = function(
  Count,
  Tcond,
  Nstep,
  Nproj,
  Tstep = 1, 
  infprofile,
  testdelayprofile,
  resultdelayalpha = 1,
  gp_time_scale = 14,
  gp_time_decay_scale = .1,
  fixed_gp_time_length_scale = -1.0,
  mu_scale = 0.5,
  sigma_scale = 1.0,
  phi_latent_scale = 5.0,
  phi_observed_scale = 5.0,
  xi_scale = 0.01,
  outlier_prob_threshold = 1.0,
  outlier_count_threshold = 10,
  reconstruct_infections = TRUE,
  num_iterations = 4000,
  num_chains = 1,
  percentiles = c(.025,.25,.5,.75,.975)
) {

  # ------------------------------------------------------------------------ #
  #                             Main computation                             #
  # ------------------------------------------------------------------------ #

  if (is.null(colnames(Count))) {
    startdate = 1
  } else if (typeof(colnames(Count)[1])=="character") {
    startdate = as.Date(colnames(Count)[1])
  } else {
    startdate = colnames(Count)[1]
  }
  Tall = length(Count)
  if (typeof(Count)=="integer") {
    dim(Count) = c(1,Tall)
  }

  data <- list(
    Count = Count,
    Tall = Tall,
    Tcond = Tcond,
    Tstep = Tstep, 
    Nstep = Nstep,
    Nproj = Nproj,
    Tip = length(infprofile),
    Ttdp = length(testdelayprofile),
    Trdp = length(resultdelayalpha),
    infprofile = infprofile,
    testdelayprofile = testdelayprofile,
    resultdelayalpha = c(resultdelayalpha,0),
    gp_time_scale = gp_time_scale,
    gp_time_decay_scale = gp_time_decay_scale,
    fixed_gp_time_length_scale = fixed_gp_time_length_scale,
    mu_scale = mu_scale,
    sigma_scale = sigma_scale,
    phi_latent_scale = phi_latent_scale,
    phi_observed_scale = phi_observed_scale,
    xi_scale = xi_scale,
    outlier_prob_threshold = outlier_prob_threshold,
    outlier_count_threshold = outlier_count_threshold,
    reconstruct_infections = reconstruct_infections
  )
  
  init = list()
  init[[1]] = list(
    mu = 0.0,
    sigma = sigma_scale,
    alpha1 = gp_time_decay_scale,
    phi_latent = phi_latent_scale,
    phi_observed = phi_observed_scale,
    xi = xi_scale
  )
  for (i in 1:Nstep) {
    init[[1]][[paste('Reta[',i,']',sep='')]] = 0.0
    init[[1]][[paste('Ceta[',i,']',sep='')]] = 0.0
  }
  notpars = c(
    "Reta",
    "Ceta",
    "Ecount"
  )
    
  
  options(mc.cores = min(num_chains,parallel::detectCores()))
  stanfit = stan(
    file = 'epimap/stan_files/Rmap-clean.stan',
    data = data, 
    init = init,
    iter = num_iterations, 
    pars = notpars,
    include = FALSE,
    chains = num_chains,
    control = list(adapt_delta = .9)
  )

  Rt = summary(stanfit, pars = "Rt", probs=percentiles)$summary
  Xt = summary(stanfit, pars = "Xt", probs=percentiles)$summary
  rownames(Rt) = as.character(seq(from=startdate+Tcond,length.out=Nstep+Nproj,by=Tstep))
  rownames(Xt) = as.character(seq(from=startdate,length.out=Tcond+Tstep*Nstep,by=1))

  list(
    stanfit = stanfit,
    Rt = Rt,
    Xt = Xt
  )
}  

