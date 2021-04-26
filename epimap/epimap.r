library(rstan)
library(gsubfn)
library(plyr)


#' Stop R without throwing errors. Allows for exiting if nothing to do in a way
#' that doesn't break bash scripts.
stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

#' Parse kernel arguments from user-friendly strings to integers for STAN.
#' 
#' @param spatialkernel The string represnetative of the spatial kernel to use. 
#' @param temporalkernel The string represnetative of the temporal kernel to use.
#' @param localkernel The string represnetative of the local kernel to use.
#' @param globalkernel The string represnetative of the global kernel to use.
#' 
#' @return A list containing the integer representations of the given kernels expected by STAN.
parse_kernels = function(
  spatialkernel,
  temporalkernel,
  localkernel,
  globalkernel
) {
  SPATIOTEMPORALKERNELS = list(
    'none' = 1,
    'expquad' = 2,
    'matern12' = 3,
    'matern32' = 4,
    'matern52' = 5
  )
  SPATIAL_KERNEL = SPATIOTEMPORALKERNELS[[opt$spatialkernel]]
  if (is.null(SPATIAL_KERNEL)) {
    stop(c('Unrecognised spatial kernel ',opt$spatialkernel));
  }
  TEMPORAL_KERNEL = SPATIOTEMPORALKERNELS[[opt$temporalkernel]]
  if (is.null(TEMPORAL_KERNEL)) {
    stop(c('Unrecognised temporal kernel ',opt$temporalkernel));
  }
  LOCALKERNELS = list(
    'none' = 0,
    'local' = 1
  )
  LOCAL_KERNEL = LOCALKERNELS[[opt$localkernel]]
  if (is.null(LOCAL_KERNEL)) {
    stop(c('Unrecognised local kernel ',opt$localkernel));
  }
  GLOBALKERNELS = list(
    'none' = 0,
    'global' = 1
  )
  GLOBAL_KERNEL = GLOBALKERNELS[[opt$globalkernel]]
  if (is.null(GLOBAL_KERNEL)) {
    stop(c('Unrecognised global kernel ',opt$globalkernel));
  }

  list(
      spatial = SPATIAL_KERNEL,
      temporal = TEMPORAL_KERNEL,
      local = LOCAL_KERNEL,
      global = GLOBAL_KERNEL
  )
}

#' Parse observation arguments from user-friendly strings to integers for STAN.
#' 
#' @param observation_data The string represnetative of the observation data to use. 
#' @param observation_model The string represnetative of the observation model to use. 
#' 
#' @return A list containing the integer representations of the given observation arguments expected by STAN.
parse_observation = function(
  observation_data,
  observation_model
) {
  OBSERVATION_DATA = list(
    'count' = 1,
    'cleaned_latent_mean' = 2,
    'cleaned_latent_sample' = 2,
    'cleaned_recon_sample' = 3,
    'latent_reports' = 4
  )
  OBSERVATION_MODEL = list(
    'poisson' = 1,
    'negative_binomial_2' = 2,
    'negative_binomial_3' = 3,
    'gaussian' = 4
  )
  OBSERVATION_DATA = OBSERVATION_DATA[[opt$observation_data]]
  if (is.null(OBSERVATION_DATA)) {
    stop(c('Unrecognised observation option ',opt$observation_data));
  }
  OBSERVATION_MODEL = OBSERVATION_MODEL[[opt$observation_model]]
  if (is.null(OBSERVATION_MODEL)) {
    stop(c('Unrecognised observation option ',opt$observation_model));
  }
  
  list(
      data = OBSERVATION_DATA,
      model = OBSERVATION_MODEL
  )
}

#' Parse metapopulation arguments from user-friendly strings to integers for STAN.
#' 
#' @param metapop "none"/"in"/"in_out".
#' @param flux A vector of matricies specifiying the fluxes between regions of differnet types.
#' @param  N The number of regoing.
#' 
#' @return A list containing the integer representations of the given observation arguments expected by STAN.
parse_metapop = function(
  metapop,
  flux,
  N
) {
  if (identical(metapop, "none")) {
      do_metapop = 0
      do_in_out = 0
  } else if (identical(metapop, "in")) {
      do_metapop = 1
      do_in_out = 0
  } else if (identical(metapop, "in_out")) {
      do_metapop = 1
      do_in_out = 1
  } else {
      stop(c('Unrecognised metapop option ',metapop));
  }

  # Fill with 
  if (is.null(flux)) {
      if (do_metapop) {
          stop(c('Doing metapop but no fluxes specified ',metapop));
      } else {
          flux = list(matrix(0.0, N, N))
      }
  }

  list(
      do = do_metapop,
      in_out = do_in_out,
      flux = flux
  )
}

#' Compute the parameter intialisation
#' 
#' 
parse_epimap_init = function(
  numchains,
  gp_space_decay,
  gp_time_decay,
  Mstep,
  Mproj,
  N,
  Tlik,
  F
) {
  lapply(1:numchains, function(i) {
    env = list2env(list(
      global_sigma = .25,
      gp_sigma = .25,
      local_scale = .25,
      local_space_sigma = .5,
      local_time_sigma = .1,
      gp_space_decay = gp_space_decay,
      gp_time_decay = gp_time_decay,
      infection_dispersion = 1.0
    ))
    setval = function(par,val,...) {
      env[[paste(par,'[',paste(...,sep=','),']',sep='')]]=val
    }
    # lapply(1:N, function(f) setval('local_space_sigma', .1, f))
    # lapply(1:Mstep+Mproj, function(f) setval('local_time_sigma', .1, f))
    lapply(1:(Mstep+Mproj), function(k) {
      lapply(1:N, function(j) {
        l = j+(k-1)*N;
        setval('local_exp', 1.0, j, k);
        lapply(
          c(
            'gp_eta_in','gp_eta_out',
            'global_eta_in','global_eta_out',
            'local_eta_in','local_eta_out'
          ),
          function(par) {
            setval(par, .1*rnorm(1,0,1) , l)
          }
        )
      })
      setval('coupling_rate',.01, k);
    })
    lapply(1:Tlik, function(s) {
      lapply(1:N, function(j) {
        setval('infection_eta', rnorm(1,0,1), s, j)
      })
    })
    lapply(1:F, function(f) setval('flux_probs', 1/F, f))
    lapply(1:N, function(j) setval('case_dispersion', 1.0, j))
    as.list(env)
  })
}

#' Run a STAN file with the given options
#' 
#' @param stan_file The STAN file to run.
#' @param stan_data The data for the STAN model.
#' @param stan_init The parameter initialisation for the STAN model.
#' @param stan_pars Parameters to output and keep.
#' @param iter The number of HMC iterations to run.
#' @param chain The number of chains to run.
#' @param control The control parameters to apply.
#' @param summary_pars The summary parameters to print.
#' 
#' @return The STAN fit object.
run_stanfile = function(
  stan_file,
  stan_data,
  stan_init,
  stan_pars,
  iter,
  chains,
  control,
  summary_pars
) {
  print(stan_file)

  options(mc.cores = min(chains,parallel::detectCores()))
  fit <- stan(
    file = stan_file,
    data = stan_data,
    init = stan_init,
    pars = stan_pars,
    iter = iter,
    chains = chains,
    control = control
  )

  print(summary(fit,
    pars=summary_pars,
    probs=c(.025,0.5,.975)
  )$summary)

  fit
}


#' Epimap inference of Rt and incidence curve.
#' 
#' Uses Bayesian inference, with a latent incidence process given by a 
#' contiuous approximation of a negative-binomial renewal process,
#' a piecwise constant prior on Rt define by a Gaussian Process, and an
#' observwaton mode for counts of positive diagnoses described by a 
#' negative-binamoial and a delay distribution. Approximates the full model
#' by only computing the model on a subset of areas, using the single area 
#' model on all other areas as an approximation for their incidence curves.
#'
#' @param iter Number of iterations of HMC to run.
#' @param chains Number of chains to run.
#' @param cases Matrix of case counts.
#' @param latent_infections Matrix of preinfered latent infections.
#' @param Tcur Index in cases/latent_infections to end conditioning on.
#' @param Tstep The length of each time step of the model.
#' @param Tpred The number of days to compute the predictive likelihood of.
#' @param Mstep The number of time steps to condition on before Tcur.
#' @param Mignore The number of time steps of data to ignore.
#' @param Mproj The number of time steps to project the model forward.
#' @param area_modelled Vector specifying if a given region should be used in modelled.
#' @param area_inferred Vector specifying if a given region should be inferred.
#' @param spatialkernel String specifying the spatial kernel.
#' @param temporalkernel String specifying the temporal kernel.
#' @param localkernel String specifying the local kernel.
#' @param globalkernel String specifying the global kernel.
#' @param observation_data String specifying the observation data to use.
#' @param observation_model String specifying the observation model to use.
#' @param infection_profile Vector specifying the infectivity profile.
#' @param delay_profile Vector specifying the observation delay profile.
#' @param geographic_distances Matrix specifiying the geographic distance between regions [areas, areas].
#' @param temporal_distances Matrix specifying the temproal distances between time steps.
#' @param temporal_correlation Matrix specifying if given timesteps should be correlated.
#' @param metapop String specifying the type of metapopulation model to use.
#' @param flux Vector of arrays containing intra regional flux. NULL / [F, [areas, areas]]
#' @param gp_space_scale
#' @param gp_space_decay_scale
#' @param gp_time_scale
#' @param gp_time_decay_scale
#' @param fixed_gp_space_length_scale If > 0.0, will use this fixed scale in place of infereing it
#' @param fixed_gp_time_length_scale If > 0.0, will use this fixed scale in place of infereing it
#' @param constant_forward_rt Whether to use the final weeks infereed Rt to project forward, or sample from the model. 1/0
#' @param full_cases_dist Whether to output the full disrtibution of projected cases, or the distribution of the mean. 1/0
#'
#' @export 
epimap_region = function(
  iter,
  chains,
  cases,
  latent_infections,
  Tcur,
  Tstep,
  Tpred,
  Mstep,
  Mignore,
  Mproj,
  area_modelled,
  area_inferred,
  spatialkernel,
  temporalkernel,
  localkernel,
  globalkernel,
  observation_data,
  observation_model,
  infection_profile,
  delay_profile,
  geographic_distances,
  temporal_distances,
  temporal_correlation,
  metapop="none",
  flux=NULL,
  gp_space_scale = 0.5,
  gp_space_decay_scale = .25,
  gp_time_scale = 50,
  gp_time_decay_scale = .25,
  fixed_gp_space_length_scale = -1.0,
  fixed_gp_time_length_scale = -1.0,
  constant_forward_rt = False,
  full_cases_dist = False
) {

  if (sum(area_inferred) == 0) {
    message("Skipping region, no areas to be inferred")
    stop_quietly()
  }

  kernels = parse_kernels(
    spatialkernel, 
    temporalkernel, 
    localkernel, 
    globalkernel
  )

  observation = parse_observation(
    observation_data,
    observation_model
  )

  metapop = parse_metapop(
    metapop, 
    flux,
    length(area_modelled)
  )
  print(cat("Mproj: ", Mproj, " Tpred: ", Tpred))
  print("Epimap data")
  epimap_data <- list(
    Nall = length(area_modelled),
    Tall = ncol(cases),
    Tcur = Tcur,
    Tstep = Tstep,
    Tpred = Tpred,
    Mstep = Mstep,
    Mignore = Mignore,
    Mproj = Mproj,

    Ct_all = cases,
    Clean_latent = latent_infections,

    area_modelled = area_modelled,
    area_inferred = area_inferred,

    Tip = length(infection_profile),
    infprofile = infection_profile,
    Tdp = length(delay_profile),
    delayprofile = delay_profile,
    F = length(metapop$flux),
    flux = metapop$flux,

    geodist_all = geographic_distances,
    timedist = temporal_distances,
    timecorcut = temporal_correlation,

    gp_space_scale = gp_space_scale,
    gp_space_decay_scale = gp_space_decay_scale,
    gp_time_scale = gp_time_scale,
    gp_time_decay_scale = gp_time_decay_scale,
    fixed_gp_space_length_scale = fixed_gp_space_length_scale,
    fixed_gp_time_length_scale = fixed_gp_time_length_scale,

    coupling_mu_loc = 0,
    coupling_mu_scale = .5,
    coupling_sigma_scale = .5,
    coupling_alpha_scale = 1-exp(-Tstep/28),
    
    SPATIAL_KERNEL = kernels$spatial,
    TEMPORAL_KERNEL = kernels$temporal,
    LOCAL_KERNEL = kernels$local,
    GLOBAL_KERNEL = kernels$global,

    DO_METAPOP = metapop$do,
    DO_IN_OUT = metapop$in_out,

    OBSERVATION_DATA = observation$data,
    OBSERVATION_MODEL = observation$model,

    CONSTANT_FORWARD_RT = constant_forward_rt,
    FULL_CASES_DISTRIBUTION = full_cases_dist
  )

  print("Epimap init")
  epimap_init = parse_epimap_init(
    chains,
    gp_space_decay_scale,
    gp_time_decay_scale,
    Mstep,
    Mproj,
    length(area_modelled),
    Mstep*Tstep,
    length(metapop$flux)
  )

  epimap_summary_pars = c(
      "global_sigma",
      "gp_sigma",
      "gp_time_length_scale",
      "gp_space_length_scale",
      "infection_dispersion",
      "case_dispersion",
      "xi",
      "infection_dispersion",
      "coupling_rate",
      "Rt_all",
      "flux_probs",
      "local_space_sigma",
      "local_time_sigma"
  )
  if (fixed_gp_space_length_scale<0.0) {
    epimap_summary_pars = c(epimap_summary_pars,"gp_space_length_scale")
  }
  if (fixed_gp_time_length_scale<0.0) {
    epimap_summary_pars = c(epimap_summary_pars,"gp_time_length_scale")
  }

  epimap_pars = c(
    epimap_summary_pars,
    "Ppred",
    "Xpred",
    "Xproj",
    "Bpred",
    "Bproj",
    "Cpred",
    "Cproj",
    "fluxproportions",
    "Rt",
    "Rt_region",
    "Bproj_region",
    "Bpred_region",
    "Cproj_region",
    "Cpred_region",
    "case_precision",
    "Xt_region",
    "Zt_region",
    "case_dispersion"
  )

  print("Epimap control")
  epimap_control = list(
    adapt_delta = .9,
    max_treedepth = 10 
  )

  fit <- run_stanfile(
    "epimap/stan_files/Rmap-latent.stan",
    epimap_data,
    epimap_init,
    epimap_pars,
    iter,
    chains,
    epimap_control,
    epimap_summary_pars
  )
  
  fit
}

#' Epimap inference of Rt and incidence curve.
#' 
#' Uses Bayesian inference, with a latent incidence process given by a 
#' contiuous approximation of a negative-binomial renewal process,
#' a piecwise constant prior on Rt define by a Gaussian Process, and an
#' observwaton mode for counts of positive diagnoses described by a 
#' negative-binamoial and a delay distribution. Approximates the full model
#' by using the incidence curves from the single area model, and inferes
#' just the Rt using the full mdoel.
#' 
#' @param iter Number of iterations of HMC to run.
#' @param chains Number of chains to run.
#' @param cases Matrix of case counts.
#' @param reconstructed_infections Matrix of preinfered reconstructed infections.
#' @param latent_infections Matrix of preinfered latent infections.
#' @param Tcur Index in cases/latent_infections to end conditioning on.
#' @param Tstep The length of each time step of the model.
#' @param Tpred The number of days to compute the predictive likelihood of.
#' @param Mstep The number of time steps to condition on before Tcur.
#' @param Mignore The number of time steps of data to ignore.
#' @param Mproj The number of time steps to project the model forward.
#' @param spatialkernel String specifying the spatial kernel.
#' @param temporalkernel String specifying the temporal kernel.
#' @param localkernel String specifying the local kernel.
#' @param globalkernel String specifying the global kernel.
#' @param observation_data String specifying the observation data to use.
#' @param observation_model String specifying the observation model to use.
#' @param sparse_region A onehot matrix specifying which areas belong to which regions. 
#' @param infection_profile Vector specifying the infectivity profile.
#' @param delay_profile Vector specifying the observation delay profile.
#' @param geographic_distances Matrix specifiying the geographic distance between regions [areas, areas].
#' @param temporal_distances Matrix specifying the temproal distances between time steps.
#' @param temporal_correlation Matrix specifying if given timesteps should be correlated.
#' @param metapop String specifying the type of metapopulation model to use.
#' @param flux Vector of arrays containing intra regional flux. NULL / [F, [areas, areas]]
#' @param gp_space_scale
#' @param gp_space_decay_scale
#' @param gp_time_scale
#' @param gp_time_decay_scale
#' @param fixed_gp_space_length_scale If > 0.0, will use this fixed scale in place of infereing it
#' @param fixed_gp_time_length_scale If > 0.0, will use this fixed scale in place of infereing it
#' @param constant_forward_rt Whether to use the final weeks infereed Rt to project forward, or sample from the model. 1/0
#' @param full_cases_dist Whether to output the full disrtibution of projected cases, or the distribution of the mean. 1/0
#'
#' @export 
epimap_twostage = function(
  iter,
  chains,
  cases,
  reconstructed_infections,
  latent_infections,
  Tcur,
  Tstep,
  Tpred,
  Mstep,
  Mignore,
  Mproj,
  spatialkernel,
  temporalkernel,
  localkernel,
  globalkernel,
  observation_data,
  observation_model,
  sparse_region,
  infection_profile,
  delay_profile,
  geographic_distances,
  temporal_distances,
  temporal_correlation,
  metapop="none",
  flux=NULL,
  gp_space_scale = 0.5,
  gp_space_decay_scale = .25,
  gp_time_scale = 50,
  gp_time_decay_scale = .25,
  fixed_gp_space_length_scale = -1.0,
  fixed_gp_time_length_scale = -1.0,
  constant_forward_rt = False,
  full_cases_dist = False
) {

  kernels = parse_kernels(
      spatialkernel, 
      temporalkernel, 
      localkernel, 
      globalkernel
  )

  observation = parse_observation(
      observation_data,
      observation_model
  )

  metapop = parse_metapop(
      metapop, 
      flux,
      length(area_modelled)
  )
  print(cat("Mproj: ", Mproj, " Tpred: ", Tpred))
  print("Epimap data")
  epimap_data <- list(
      Nall = dim(sparse_region)[1],
      Tall = ncol(cases),
      Tcur = Tcur,
      Tstep = Tstep,
      Tpred = Tpred,
      Mstep = Mstep,
      Mignore = Mignore,
      Mproj = Mproj,

      Ct_all = cases,
      Clean_recon = reconstructed_infections,
      Clean_latent = latent_infections,

      Tip = length(infection_profile),
      infprofile = infection_profile,
      Tdp = length(delay_profile),
      delayprofile = delay_profile,
      F = length(metapop$flux),
      flux = metapop$flux,

      N_region = dim(sparse_region)[2],
      sparse_region = sparse_region,

      geodist_all = geographic_distances,
      timedist = temporal_distances,
      timecorcut = temporal_correlation,

      gp_space_scale = gp_space_scale,
      gp_space_decay_scale = gp_space_decay_scale,
      gp_time_scale = gp_time_scale,
      gp_time_decay_scale = gp_time_decay_scale,
      fixed_gp_space_length_scale = fixed_gp_space_length_scale,
      fixed_gp_time_length_scale = fixed_gp_time_length_scale,

      coupling_mu_loc = -2.19, # centre at .1
      coupling_mu_scale = 0.25, # set mean of process to be 0.1, 1 std = 0.024-0.33
      coupling_sigma_scale = .25,
      coupling_alpha_scale = 1-exp(-Tstep/28),
      
      SPATIAL_KERNEL = kernels$spatial,
      TEMPORAL_KERNEL = kernels$temporal,
      LOCAL_KERNEL = kernels$local,
      GLOBAL_KERNEL = kernels$global,

      DO_METAPOP = metapop$do,
      DO_IN_OUT = metapop$in_out,

      OBSERVATION_DATA = observation$data,
      OBSERVATION_MODEL = observation$model,

      CONSTANT_FORWARD_RT = constant_forward_rt,
      FULL_CASES_DISTRIBUTION = full_cases_dist
  )

  print("Epimap init")
  epimap_init = parse_epimap_init(
      chains,
      gp_space_decay_scale,
      gp_time_decay_scale,
      Mstep,
      Mproj,
      dim(sparse_region)[1],
      Mstep*Tstep,
      length(metapop$flux)
  )

  epimap_summary_pars = c(
      "global_sigma",
      "gp_sigma",
      "gp_time_length_scale",
      "gp_space_length_scale",
      # "phi_latent",
      # "phi_observed",
      "xi",
      "infection_dispersion",
      "coupling_rate",
      "Rt_all",
      "flux_probs"
  )
  if (fixed_gp_space_length_scale<0.0) {
    epimap_summary_pars = c(epimap_summary_pars,"gp_space_length_scale")
  }
  if (fixed_gp_time_length_scale<0.0) {
    epimap_summary_pars = c(epimap_summary_pars,"gp_time_length_scale")
  }

  epimap_pars = c(
      epimap_summary_pars,
      "Ppred",
      "Bpred",
      "Bproj",
      "Cpred",
      "Cproj",
      "Xpred",
      "Xproj",
      "fluxproportions",
      "Rt",
      "Rt_region",
      "Cproj_region",
      "Cpred_region",
      "Bproj_region",
      "Bpred_region"
  )

  print("Epimap control")
  epimap_control = list(
      adapt_delta = .9,
      max_treedepth = 10 
  )

  fit <- run_stanfile(
      "epimap/stan_files/Rmap.stan",
      epimap_data,
      epimap_init,
      epimap_pars,
      iter,
      chains,
      epimap_control,
      epimap_summary_pars
  )
  
  fit
}



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
#' @param iter Number of MCMC iterations; default 3000
#' @param chains Number of MCMC chains; default 1
#' @param cases Dataframe or list of case counts of positive tests by day; colnames(Count) can be date strings in YYYY-MM-DD format.
#' @param Tcond Initial number of days that are not modelled (but are conditioned on)
#' @param Nstep Number of AR1 time steps to model 
#' @param Nproj Number of AR1 time steps to forecast
#' @param Tpred Number of days after the period used for the liklehood to compute the predictive likelihood for 
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
#' @param phi_observed_scale Scale of prior on dispersion parameter of observation process of positive test counts; default 10.0
#' @param xi_scale Scale of prior on exegeneous infection rate
#' @param reconstruct_infections Whether to reconstruct infection day for each case; default True
#' @param outlier_prob_threshold Probability hreshold to determine if a diagnosis count is an outlier; default 1.0
#' @param outlier_count_threshold Count threshold to determine if a diagnosis count is an outlier; default 10
#'
#' @return A list consists of the stanfit MCMC output.
#'
#' @export
epimap_singlearea = function(
  iter,
  chains,
  cases,
  Tcond,
  Nstep,
  Nproj,
  Tpred,
  Tstep = 1,
  infprofile,
  testdelayprofile,
  resultdelayalpha = 1,
  gp_time_scale = 14,
  gp_time_decay_scale = .1,
  fixed_gp_time_length_scale = -1.0,
  mu_scale = 0.5,
  sigma_scale = 1.0,
  phi_latent_scale = 10.0,
  phi_observed_scale = 10.0,
  xi_scale = 0.01,
  outlier_prob_threshold = 1.0,
  outlier_count_threshold = 10,
  reconstruct_infections = TRUE
) {

  # ------------------------------------------------------------------------ #
  #                             Main computation                             #
  # ------------------------------------------------------------------------ #

  if (is.null(colnames(cases))) {
    startdate = 1
  } else if (typeof(colnames(cases)[1])=="character") {
    startdate = as.Date(colnames(cases)[1])
  } else {
    startdate = colnames(cases)[1]
  }
  Tall = length(cases)
  if (typeof(cases)=="integer") {
    dim(cases) = c(1,Tall)
  }

  print(cat("Nproj: ", Nproj, " Tpred: ", Tpred))
  print("Epimap data")
  epimap_data <- list(
    Count = cases,
    Tall = Tall,
    Tcond = Tcond,
    Tstep = Tstep, 
    Nstep = Nstep,
    Nproj = Nproj,
    Tpred = Tpred,
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
  
  print("Epimap init")
  epimap_init = list()
  epimap_init[[1]] = list(
    mu = 0.0,
    sigma = sigma_scale,
    alpha1 = gp_time_decay_scale,
    phi_latent = phi_latent_scale,
    phi_observed = phi_observed_scale,
    xi = xi_scale
  )
  for (i in 1:Nstep) {
    epimap_init[[1]][[paste('Reta[',i,']',sep='')]] = 0.0
    epimap_init[[1]][[paste('Ceta[',i,']',sep='')]] = 0.0
  }

  epimap_summary_pars = c(
    "mu",
    "sigma",
    "alpha",
    "gp_time_length_scale",
    "phi_latent",
    "phi_observed",
    "xi",
    "Noutliers",
    "meandelay",
    "resultdelayprofile",
    "Rx",
    "Rt"
  )

  epimap_pars = c(
    epimap_summary_pars,
    "Crecon",
    "Cpred",
    "Cproj",
    "Bpred",
    "Bproj",
    "Xpred",
    "Xproj",
    "Noutliers",
    "Xt",
    "Xt_proj",
    "Ppred"
  )

  print("Epimap control")
  epimap_control = list(
      adapt_delta = .9,
      max_treedepth = 10 
  )

  fit <- run_stanfile(
    'epimap/stan_files/Rmap-clean.stan',
    epimap_data,
    epimap_init,
    epimap_pars,
    iter,
    chains,
    epimap_control,
    epimap_summary_pars
  )

  fit
}  
