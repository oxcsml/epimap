library(rstan)
library(geosphere)
library(optparse)
source('dataprocessing/read_data.r')
source('covidmap/utils.r')

Rmap_options = function(
  spatialkernel        = "none",
  temporalkernel       = "matern12",
  localkernel          = "local",
  globalkernel         = "none",
  gp_space_scale       = 0.5, # units of 100km
  gp_space_decay_scale = .25,
  gp_time_scale        = 50.0, # units of 1 day
  gp_time_decay_scale  = .25,
  fixed_gp_space_length_scale = -1, # 0.5,
  fixed_gp_time_length_scale = -1, # 100.0,
  constant_forward_rt  = 0, # if 0, use the Rt from last week to predict forwards
  full_cases_distribution = 1,
  metapop              = "traffic_forward,traffic_reverse,uniform,in",
  #metapop              = "traffic_forward,traffic_reverse,radiation1,radiation2,radiation3,uniform,in",
  observation_data     = "cleaned_latent_sample",
  observation_model    = "gaussian",
  cleaned_sample_id    = 1, 

  first_day_modelled   = NULL,
  last_day_modelled    = NULL,
  weeks_modelled       = 15,
  days_ignored         = 7,
  days_per_step        = 7,
  days_predicted       = 2,
  steps_ignored_stage2 = 1,
  num_steps_forecasted = 3,

  thinning             = 10,
  chains               = 1,
  iterations           = 3000, 

  data_directory       = "data/",
  clean_directory      = "fits/clean", 
  results_directory    = NULL,

  limit_area           = NULL,
  limit_radius         = NULL,

  stage                = "map"
) {
  as.list(environment())
}

##########################################################################
##########################################################################
Rmap_setup = function(opt = Rmap_options()) {
  env = new.env(parent=globalenv())
  env$opt = opt
  Rmap_read_data(env)
  with(env,{

    #########################################################
    numchains = opt$chains
    numiters = opt$iterations

    options(mc.cores = min(numchains,parallel::detectCores()))
    rstan_options(auto_write = TRUE)

    #########################################################
    list[Mstep, Tstep, Tcond, Tlik, Tcur, Tignore] = process_dates_modelled(
      dates,
      first_day_modelled = opt$first_day_modelled,
      last_day_modelled  = opt$last_day_modelled,
      days_ignored       = opt$days_ignored,
      weeks_modelled     = opt$weeks_modelled,
      days_per_step      = opt$days_per_step
    )
    Tpred = opt$days_predicted    # number of days for predictive probs eval
    Mproj = opt$num_steps_forecasted
    Tproj = Tstep*Mproj           # number of days to project forward

    message("Tcur: ", Tcur, ", length Clean_latent: ", length(Clean_latent))
    stopifnot(Tcur == length(Clean_latent))
      message("Tcur: ", Tcur, ", length Clean_recon: ", length(Clean_recon))
    stopifnot(Tcur == length(Clean_recon))

    days_likelihood = dates[(Tcond+1):Tcur]
    days_pred_held_out = seq(dates[Tcur+1],by=1,length.out=Tpred)
    days_proj = seq(dates[Tcur+1],by=1,length.out=Tproj)
    days_forward = seq(dates[Tcur+1],by=1,length.out=max(Tproj, Tpred))
    days = c(days_likelihood, days_forward)
    # days = c(days_likelihood)
    # Mproj = 0


    message("Days used for held out likelihood = ",
      days_pred_held_out[1],"...",days_pred_held_out[Tpred])

    if (is.null(opt$results_directory)) {
      opt$results_directory = paste(
        'fits/',
        as.character(Sys.time(),format='%Y%m%d'),
        '-',as.character(days_likelihood[length(days_likelihood)],format='%Y%m%d'),
        '-',opt$spatialkernel,
        '-',opt$temporalkernel,
        '-',opt$localkernel,
        '-',opt$globalkernel,
        '-',opt$metapop,
        '-',opt$observation_data,
        '-',opt$observation_model,
        sep=''
      )
    }
    dir.create(opt$results_directory, showWarnings = FALSE)
    if (opt$cleaned_sample_id>0) {
      runname = paste(opt$results_directory,'/',opt$cleaned_sample_id,'_',sep='')
    } else {
      runname = paste(opt$results_directory,'/',sep='')
    }
    message("runname = ",runname)

    #########################################################
  })
  env # return environment
} # Rmap_setup
##########################################################################
##########################################################################

##########################################################################
##########################################################################
Rmap_run = function(env) {
  with(env, {

    start_time <- Sys.time()

    #########################################################
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


    #########################################################
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

    #########################################################
    # metapopulation cross-area fluxes.
    METAPOPMODEL = strsplit(opt$metapop,',')[[1]]
    METAPOPOPTIONS = list(
      'radiation1' = radiation_flux[,,1], # smoothed radiation model with length scale = .1 (10km)
      'radiation2' = radiation_flux[,,2], # smoothed radiation model with length scale = .2 (20km)
      'radiation3' = radiation_flux[,,3], # smoothed radiation model with length scale = .5 (50km)
      'traffic_forward' = traffic_flux[,,1], # infected commuters taking infection from home to work
      'traffic_reverse' = traffic_flux[,,2], # commuters getting infected at work and bringing back home
      'uniform' = matrix(1.0/N,N,N) # uniform cross-area infection
    )
    if (length(METAPOPMODEL)==1 && METAPOPMODEL[1] == 'none') {
      DO_METAPOP = 0
      DO_IN_OUT = 0
      flux = array(0,dim=c(0,N,N));
      F = 0
    } else {
      DO_METAPOP = 1
      if ('in_out' %in% METAPOPMODEL) {
        DO_IN_OUT = 1
      } else if ('in' %in% METAPOPMODEL) {
        DO_IN_OUT = 0
      } else {
        stop(c('One of "in" or "in_out" should be in metapop option ',opt$metapop));
      }
      flux = list()
      F = 0
      for (s in METAPOPMODEL) {
        if (! s %in% c('in','in_out')) {
          t = METAPOPOPTIONS[[s]]
          if (is.null(t)) {
            stop(c('Unrecognised metapop option in ',opt$observation));
          }
          F = F+1
          flux[[F]] = t
        }
      }
    }

    #########################################################
    # time steps
    times = 1:(Mstep+Mproj)
    timedist = matrix(0, Mstep+Mproj, Mstep+Mproj)
    for (i in 1:(Mstep+Mproj)) {
      for (j in 1:(Mstep+Mproj)) {
        timedist[i, j] = abs(times[i] - times[j]) * Tstep
      }
    }

    # precompute lockdown cutoff kernel
    lockdown_day = as.Date("2020-03-23")
    days_period_start = days[seq(1, length(days), Tstep)]
    days_period_start = vapply(days_period_start, (function (day) as.Date(day, format="%Y-%m-%d")), double(1))
    day_pre_lockdown = vapply(days_period_start, (function (day) day < lockdown_day), logical(1))
    
    time_corellation_cutoff = matrix(0,Mstep+Mproj,Mstep+Mproj)
    for (i in 1:(Mstep+Mproj)) {
      for (j in 1:(Mstep+Mproj)) {
        time_corellation_cutoff[i, j] = !xor(day_pre_lockdown[i], day_pre_lockdown[j])
      }
    }

    #########################################################
    # Main computation
    print(cat("Mproj: ", Mproj, " Tpred: ", Tpred))
    Rmap_data <- list(
      N = N,
      Mstep = Mstep,
      Tall = Tall,
      Tcur = Tcur,
      Tstep = Tstep,
      Mignore = opt$steps_ignored_stage2,
      Mproj = Mproj,
      Tpred = Tpred,

      Ct = AllCount,
      Clean_latent = Clean_latent,
      Clean_recon = Clean_recon,
      geodist = geodist,
      timedist = timedist,
      timecorcut = time_corellation_cutoff,

      gp_space_scale = opt$gp_space_scale,
      gp_space_decay_scale = opt$gp_space_decay_scale,
      gp_time_scale = opt$gp_time_scale,
      gp_time_decay_scale = opt$gp_time_decay_scale,
      fixed_gp_space_length_scale = opt$fixed_gp_space_length_scale,
      fixed_gp_time_length_scale = opt$fixed_gp_time_length_scale,
      coupling_mu_loc = -2.19, # centre at .1
      coupling_mu_scale = 0.25, # set mean of process to be 0.1, 1 std = 0.024-0.33
      coupling_sigma_scale = 0.1,
      coupling_alpha_scale = 1-exp(-Tstep/28),
      
      SPATIAL_KERNEL = SPATIAL_KERNEL,
      TEMPORAL_KERNEL = TEMPORAL_KERNEL,
      LOCAL_KERNEL = LOCAL_KERNEL,
      GLOBAL_KERNEL = GLOBAL_KERNEL,
      DO_METAPOP = DO_METAPOP,
      DO_IN_OUT = DO_IN_OUT,
      OBSERVATION_DATA = OBSERVATION_DATA,
      OBSERVATION_MODEL = OBSERVATION_MODEL,
      CONSTANT_FORWARD_RT = opt$constant_forward_rt,
      FULL_CASES_DISTRIBUTION = opt$full_cases_distribution,

      Tip = Tip,
      infprofile = infprofile,
      Tdp = Tdp,
      delayprofile = testdelayprofile,
      F = F,
      flux = flux,

      N_region = N_region,
      sparse_region = sparse_region
    )

    #########################################################
    Rmap_init = lapply(1:numchains, function(i) {
      env = list2env(list(
        global_sigma = .25,
        gp_sigma = .25,
        local_space_sigma = .25,
        local_time_sigma = .1,
        gp_space_decay = opt$gp_space_decay,
        gp_time_decay = opt$gp_time_decay,
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
              # 'global_eta_in','global_eta_out',
              'local_eta_in','local_eta_out'
            ),
            function(par) {
              setval(par, rnorm(1,0,1) , l)
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
      lapply(1:N, function(j) setval('case_precision', 0.1, j))
      as.list(env)
    })
    # print(Rmap_init)
    # message(Rmap_init[[1]])
    Rmap_summary_pars = c(
      "global_sigma",
      "local_space_sigma",
      "local_time_sigma",
      "gp_sigma",
      "infection_dispersion",
      "coupling_rate",
      "flux_probs",
      "Rt",
      "Rt_all"
    )
    if (opt$fixed_gp_space_length_scale<0.0) {
      Rmap_summary_pars = c(Rmap_summary_pars,"gp_space_length_scale")
    }
    if (opt$fixed_gp_time_length_scale<0.0) {
      Rmap_summary_pars = c(Rmap_summary_pars,"gp_time_length_scale")
    }
    Rmap_pars = c(Rmap_summary_pars,
      "case_precision",
      "Ppred",
      "Cpred",
      "Cproj",
      "fluxproportions"
      # "Rt_region",
      # "Cproj_region",
      # "Cpred_region"
    )
    Rmap_control = list(
      # max_treedepth = 3, # testing only
      adapt_delta = .9,
      max_treedepth = 10 
    )

    #########################################################
    fit <- stan(
      file = 'mapping/stan_files/Rmap-latent.stan',
      data = Rmap_data,
      init = Rmap_init,
      pars = Rmap_pars,
      iter = numiters,
      chains = numchains,
      control = list(adapt_delta = .9)
    )

    #########################################################
    saveRDS(fit, paste(runname,'stanfit.rds', sep=''))

    #########################################################
    # Summary of fit
    print(summary(fit,
        pars=Rmap_summary_pars,
        probs=c(.025,0.5,.975)
    )$summary)


    #########################################################
    end_time <- Sys.time()
    message(paste(
      "Time to run:",
      round(difftime(end_time, start_time, units="hours"),3),
      "hours")
    )
    message(runname)

    #########################################################
  })
  env
} # Rmap_run
##########################################################################
##########################################################################


##########################################################################
##########################################################################
Rmap_load = function(env) {
  with(env, {
    fit = readRDS(paste(runname,'stanfit.rds',sep=''))
  })
  env
} # Rmap_load
##########################################################################
##########################################################################

##########################################################################
##########################################################################
Rmap_postprocess = function(env) {
  with(env, {
    writeresults = function(data,filename,...) {
      write.csv(data,sprintf('%s%s.csv',runname,filename),...)
    }

    # save raw samples
    save_samples = function(pars,N_sites=N,areafirst=FALSE) {
      samples = extract(fit,pars=pars,permuted=FALSE)
      samples = samples[seq(numchains,by=numchains,to=numiters/2),,,drop=FALSE]
      ns = numiters/2
      np = dim(samples)[3]/N_sites
      parnames = dimnames(samples)[[3]]
      dim(samples) = c(ns,np*N_sites)
      colnames(samples) = parnames
      if (areafirst) {
        ind = N_sites*(rep(1:np,N_sites)-1) + rep(1:N_sites,each=np)
        samples = samples[,ind,drop=FALSE]
      } else {
        samples = samples[,,drop=FALSE]
      }
      samples = round(samples,2)
      writeresults(samples,paste(pars,'_samples',sep=''))
    }
    save_samples("Rt")
    save_samples("Cpred",areafirst=TRUE)
    save_samples("Cproj",areafirst=TRUE)
    #save_samples("Rt_region", N_sites=N_region)
    #save_samples("Cpred_region", N_sites=N_region, areafirst=TRUE)
    #save_samples("Cproj_region", N_sites=N_region, areafirst=TRUE)
    #################################################################
    area_date_dataframe <- function(areas,dates,provenance,data,data_names) {
      numareas <- length(areas)
      numdates <- length(dates)
      dates <- rep(dates,numareas)
      dim(dates) <- c(numareas*numdates)
      provenance <- rep(provenance,numareas)
      dim(provenance) <- c(numareas*numdates)
      areas <- rep(areas,numdates)
      dim(areas) <- c(numareas,numdates)
      areas <- t(areas)
      dim(areas) <- c(numareas*numdates)
      df <- data.frame(area=areas,Date=dates,data=data,provenance=provenance)
      colnames(df)[3:(ncol(df)-1)] <- data_names
      df
    }

    provenance <- c(rep('inferred',Tlik),rep('projected',Tproj))
    days_all <- c(days_likelihood,seq(days_likelihood[Tlik]+1,by=1,length.out=Tproj))

    #################################################################
    # Rt posterior
    s <- summary(fit, pars="Rt", probs=c(0.025, .1, .2, 0.25, .3, .4, .5, .6, .7, 0.75, .8, .9, .975))$summary
    Rt <- s[,c("2.5%","10%", "20%", "25%", "30%", "40%", "50%", "60%", "70%", "75%", "80%", "90%","97.5%")]
    #s <- summary(fit, pars="Rt", probs=c(.1, .2, .3, .4, .5, .6, .7, .8, .9))$summary
    #Rt <- s[,c("10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%")]
    
    times = rep(c(1:(Mstep+Mproj)), N)
    places = rep(1:N, each=Mstep+Mproj)
    indicies = places + (N)*(times-1)
    Rt = Rt[indicies,]

    Rt <- Rt[sapply(1:(N*(Mstep+Mproj)),function(i)rep(i,Tstep)),]
    message(sprintf("median Rt range: [%f, %f]",min(Rt[,"50%"]),max(Rt[,"50%"])))
    df <- area_date_dataframe(
      quoted_areas,
      days_all,
      provenance,
      format(round(Rt,2),nsmall=2),
      #c("Rt_10","Rt_20","Rt_30","Rt_40","Rt_50","Rt_60","Rt_70","Rt_80","Rt_90")
      c("Rt_2_5","Rt_10","Rt_20","Rt_25","Rt_30","Rt_40","Rt_50",
        "Rt_60","Rt_70","Rt_75","Rt_80","Rt_90","Rt_97_5")
    )
    writeresults(df, 'Rt', row.names=FALSE, quote=FALSE)


    #################################################################
    # Rt exceedance probabilities
    thresholds = c(.8, .9, 1.0, 1.1, 1.2, 1.5, 2.0)
    numthresholds = length(thresholds)
    numsamples = numchains*numiters/2
    Rt <- as.matrix(fit, pars="Rt")
    dim(Rt) <- c(numsamples,Mstep+Mproj,N)
    Pexceedance = array(0.0,dim=c(Mstep+Mproj,N,numthresholds))
    for (k in 1:(Mstep+Mproj)) {
      for (i in 1:N) {
        for (x in 1:numthresholds) {
          Pexceedance[k,i,x] = mean(Rt[,k,i]>thresholds[x])
        }
      }
    }
    Pexceedance = Pexceedance[c(1:(Mstep+Mproj)),,]
    Pexceedance <- Pexceedance[sapply(1:(Mstep+Mproj),function(k)rep(k,Tstep)),,]
    dim(Pexceedance) <- c(Tstep*(Mstep+Mproj)*N,numthresholds)
    df <- area_date_dataframe(
        quoted_areas,
        days_all,
        provenance,
        format(round(Pexceedance,2),nsmall=2),
        c("P_08","P_09","P_10","P_11","P_12","P_15","P_20")
    )
    writeresults(df, 'Pexceed', row.names=FALSE, quote=FALSE)




    #################################################################
    # posterior predictives and projections
    s <- summary(fit, pars=c("Cpred"), probs=c(0.025, .25, .5, .75, .975))$summary
    Cpred <- s[,c("2.5%","25%", "50%","75%", "97.5%")]
    #Cpred <- s[,c("2.5%", "50%", "97.5%")]
    Cpred <- t(t(Cpred))
    message(sprintf("median Cpred range: [%f, %f]",min(Cpred[,"50%"]),max(Cpred[,"50%"])))
    df <- area_date_dataframe(
      quoted_areas,
      seq(dates[Tcond]+1,by=1,length.out=Mstep*Tstep),
      rep('inferred',Tlik),
      format(round(Cpred,1),nsmall=1),
      #c("C_2_5","C_25","C_50","C_75","C_97_5")
      c("C_025","C_25","C_50","C_75","C_975")
    )
    writeresults(df, 'Cpred', row.names=FALSE, quote=FALSE)


    # weekly counts. Includes 1 last column of actual counts among days ignored in model
    Tweek = Tstep # assumes Tstep = 7
    Cweekly <- as.matrix(AllCount[,(Tcond+1):(Tcond+Tlik)])
    dim(Cweekly) <- c(N,Tstep,Mstep)
    Cweekly <- apply(Cweekly,c(1,3),sum)

    # ignoredweek <- apply(AllCount[,(Tcur+1):(Tcur+Tstep)],c(1),sum)
    # Cweekly <- cbind(Cweekly,ignoredweek)

    s <- summary(fit, pars=c("Cproj"), probs=c(.5))$summary
    projectedweeks <- as.matrix(s[,"50%"])
    dim(projectedweeks) <- c(Tstep,Mproj,N)
    projectedweeks <- projectedweeks[,1:Mproj,,drop=FALSE]
    projectedweeks <- apply(projectedweeks,c(2,3),sum)
    projectedweeks <- t(projectedweeks)

    Cweekly <- cbind(Cweekly,projectedweeks)
    Cweekly <- t(Cweekly)
    Cweekly <- Cweekly[sapply(1:(Mstep+Mproj),function(k)rep(k,Tstep)),]
    dim(Cweekly) <- c(N*(Tlik+Tstep*Mproj))

    Cweekly_provenance <- c(rep('actual',Tlik),rep('projected',Tproj))
    df <- area_date_dataframe(
      quoted_areas,
      days_all,
      provenance,
      format(round(Cweekly,1),digits=6),
      c("C_weekly")
    )
    writeresults(df, 'Cweekly', row.names=FALSE, quote=FALSE)



    s <- summary(fit, pars=c("Cproj"), probs=c(0.025, .25, .5, .75, .975))$summary
    Cproj <- s[,c("2.5%","25%", "50%","75%", "97.5%")]
    #Cproj <- s[,c("2.5%", "50%", "97.5%")]
    Cproj <- t(t(Cproj))
    message(sprintf("median Cproj range: [%f, %f]",min(Cproj[,"50%"]),max(Cproj[,"50%"])))
    df <- area_date_dataframe(
      quoted_areas,
      seq(dates[Tcur]+1,by=1,length.out=Tproj),
      rep('projected',Tproj),
      format(round(Cproj,1),digits=5),
      #c("C_2_5","C_25","C_50","C_75","C_97_5")
      c("C_025","C_25","C_50","C_75","C_975")
    )
    writeresults(df, 'Cproj', row.names=FALSE, quote=FALSE)


    #################################################################
    # predictive probabilities
    s <- summary(fit, pars="Ppred", probs=c(0.025, .5, .975))$summary
    Ppred <- s[,"mean"]
    logpred <- log(Ppred)
    dim(logpred) <- c(Tpred,N)
    logpred <- t(logpred)
    message(sprintf("mean log predictives = %f",mean(logpred)))
    df <- data.frame(area = quoted_areas, logpred = logpred, provenance=rep('inferred', length(quoted_areas)))
    for (i in 1:Tpred)
      colnames(df)[i+1] <- sprintf('logpred_day%d',i)
    writeresults(df, 'logpred', row.names=FALSE, quote=FALSE)


    ####################################################################
    # pairs plot
    #pdf(paste('fits/',runname,'_pairs.pdf',sep=''),width=9,height=9)
    #pairs(fit, pars=c(
    #    "gp_space_length_scale","gp_sigma","gp_time_length_scale",
    #    "global_sigma","local_sigma","precision"))
    #dev.off()

    ####################################################################
  })
  env
} # Rmap_postprocess
##########################################################################
##########################################################################

Rmap_merge = function(env,cleaned_sample_ids) {
env$cleaned_sample_ids = cleaned_sample_ids
with(env, {

writemergedresults = function(data,filename,...) {
  write.csv(data,sprintf('%s/merged_%s.csv',opt$results_directory,filename),...)
}

numruns = length(cleaned_sample_ids)
load_samples = function(pars) {
  samples = do.call(rbind, lapply(1:numruns, function(i) {
    read.csv(paste(
      opt$results_directory,
      '/',cleaned_sample_ids[i],
      '_',pars,
      '_samples.csv',
      sep=''
    ))
  }))
  samples[,2:ncol(samples)]
}


#################################################################

area_date_dataframe <- function(areas,dates,provenance,data,data_names) {
  numareas <- length(areas)
  numdates <- length(dates)
  dates <- rep(dates,numareas)
  dim(dates) <- c(numareas*numdates)
  provenance <- rep(provenance,numareas)
  dim(provenance) <- c(numareas*numdates)
  areas <- rep(areas,numdates)
  dim(areas) <- c(numareas,numdates)
  areas <- t(areas)
  dim(areas) <- c(numareas*numdates)
  df <- data.frame(area=areas,Date=dates,data=data,provenance=provenance)
  colnames(df)[3:(ncol(df)-1)] <- data_names
  df
}


provenance <- c(rep('inferred',Tlik),rep('projected',Tproj))
days_all <- c(days_likelihood,seq(days_likelihood[Tlik]+1,by=1,length.out=Tproj))

#################################################################
# Rt posterior
Rt_samples = load_samples('Rt')
# TODO: we get very infrequent Nas in the cori model
# if (any(is.na(Rt_samples))) {
#   message("WARNING: NAs in Rt samples")
#   Rt_samples = Rt_samples[complete.cases(Rt_samples),]
# }

Rt = t(apply(Rt_samples,2,quantile,
    probs=c(0.025, .1, .2, 0.25, .3, .4, .5, .6, .7, 0.75, .8, .9, .975)
))

Rt = Rt[sapply(1:N,function(i)rep((i-1)*(Mstep+Mproj)+c(1:(Mstep+Mproj)),each=Tstep)),]
#Rt = Rt[sapply(1:N,function(i)rep((i-1)*Mstep+c(1:Mstep,rep(Mstep,Mproj)),each=Tstep)),]
df <- area_date_dataframe(
    quoted_areas,
    days_all,
    provenance,
    format(round(Rt,2),nsmall=2),
    #c("Rt_10","Rt_20","Rt_30","Rt_40","Rt_50","Rt_60","Rt_70","Rt_80","Rt_90")
    c("Rt_2_5","Rt_10","Rt_20","Rt_25","Rt_30","Rt_40","Rt_50",
      "Rt_60","Rt_70","Rt_75","Rt_80","Rt_90","Rt_97_5")
)
writemergedresults(df, 'Rt', row.names=FALSE, quote=FALSE)
message('done Rt')

#################################################################
# Rt exceedance probabilities
thresholds = c(.8, .9, 1.0, 1.1, 1.2, 1.5, 2.0)
numthresholds = length(thresholds)
numsamples = numruns * floor(numiters/2 / opt$thinning)
# numsamples = dim(env$Rt_samples)[1] 
Rt <- as.matrix(Rt_samples)
dim(Rt) <- c(numsamples,Mstep+Mproj,N)
Pexceedance = array(0.0,dim=c(Mstep+Mproj,N,numthresholds))
for (k in 1:(Mstep+Mproj)) {
  for (i in 1:N) {
    for (x in 1:numthresholds) {
      Pexceedance[k,i,x] = mean(Rt[,k,i]>thresholds[x])
    }
  }
}
Pexceedance = Pexceedance[c(1:(Mstep+Mproj)),,]
Pexceedance <- Pexceedance[sapply(1:(Mstep+Mproj),function(k)rep(k,Tstep)),,]
dim(Pexceedance) <- c(Tstep*(Mstep+Mproj)*N,numthresholds)
df <- area_date_dataframe(
    quoted_areas,
    days_all,
    provenance,
    format(round(Pexceedance,2),nsmall=2),
    c("P_08","P_09","P_10","P_11","P_12","P_15","P_20")
)
writemergedresults(df, 'Pexceed', row.names=FALSE, quote=FALSE)
message('done Pexceedance')

rm(Rt_samples)

#################################################################
# posterior predictives and projections
Cpred_samples = load_samples('Cpred')
Cpred = t(apply(Cpred_samples,2,quantile,
    probs=c(0.025, 0.25, .5, 0.75, .975)
))

df <- area_date_dataframe(
    quoted_areas,
    seq(dates[Tcond]+1,by=1,length.out=Mstep*Tstep),
    rep('inferred',Tlik),
    format(round(Cpred,1),nsmall=1),
    #c("C_2_5","C_25","C_50","C_75","C_97_5")
    c("C_025","C_25","C_50","C_75","C_975")
)
writemergedresults(df, 'Cpred', row.names=FALSE, quote=FALSE)
message('done Cpred')

rm(Cpred_samples)

####################################################################################
# weekly counts. Includes 1 last column of actual counts among days ignored in model
Cproj_samples = load_samples('Cproj')
Cweekly <- as.matrix(AllCount[,(Tcond+1):(Tcond+Tlik)])
dim(Cweekly) <- c(N,Tstep,Mstep)
Cweekly <- apply(Cweekly,c(1,3),sum)

stopifnot(Tcur+Tstep<=length(AllCount))

ignoredweek <- apply(AllCount[,(Tcur+1):(Tcur+Tstep)],c(1),sum)
Cweekly <- cbind(Cweekly,ignoredweek)

projectedweeks = as.matrix(apply(Cproj_samples,2,quantile,
    probs=c(.5)
))
dim(projectedweeks) <- c(Tstep,Mproj,N)
projectedweeks <- projectedweeks[,1:Mproj,,drop=FALSE]
projectedweeks <- apply(projectedweeks,c(2,3),sum)
projectedweeks <- t(projectedweeks)
Cweekly <- cbind(Cweekly,projectedweeks)

Cweekly <- t(Cweekly)
Cweekly <- Cweekly[sapply(1:(Mstep+Mproj),function(k)rep(k,Tstep)),]
dim(Cweekly) <- c(N*(Tlik+Tstep*Mproj))
df <- area_date_dataframe(
    quoted_areas,
    days_all,
    provenance,
    format(Cweekly,digits=3),
    c("C_weekly")
)
writemergedresults(df, 'Cweekly', row.names=FALSE, quote=FALSE)
message('done Cweekly')



Cproj = t(apply(Cproj_samples,2,quantile,
    probs=c(.025,.25,.5,.75,.975)
))
df <- area_date_dataframe(
    quoted_areas,
    seq(dates[Tcur]+1,by=1,length.out=Tproj),
    rep('projected',Tproj),
    format(round(Cproj,1),nsmall=1),
    #c("C_2_5","C_25","C_50","C_75","C_97_5")
    c("C_025","C_25","C_50","C_75","C_975")
)
writemergedresults(df, 'Cproj', row.names=FALSE, quote=FALSE)
message('done Cproj')
# Rt_region posterior
Rt_region_samples = load_samples('Rt_region')
Rt_region = t(apply(Rt_region_samples,2,quantile,
    probs=c(0.025, .1, .2, 0.25, .3, .4, .5, .6, .7, 0.75, .8, .9, .975)
))

Rt_region = Rt_region[sapply(1:N_region,function(i)rep((i-1)*Mstep+c(1:Mstep,rep(Mstep,Mproj)),each=Tstep)),]
df <- area_date_dataframe(
    quoted_regions,
    days_all,
    provenance,
    format(round(Rt_region,2),nsmall=2),
    #c("Rt_10","Rt_20","Rt_30","Rt_40","Rt_50","Rt_60","Rt_70","Rt_80","Rt_90")
    c("Rt_2_5","Rt_10","Rt_20","Rt_25","Rt_30","Rt_40","Rt_50",
      "Rt_60","Rt_70","Rt_75","Rt_80","Rt_90","Rt_97_5")
)
writemergedresults(df, 'Rt_region', row.names=FALSE, quote=FALSE)
message('done Rt_region')
# Cproj_region posterior predictive
Cproj_region_samples = load_samples('Cproj_region')
Cproj_region = t(apply(Cproj_region_samples,2,quantile,
    probs=c(.025,.25,.5,.75,.975)
))
df <- area_date_dataframe(
    quoted_regions,
    seq(dates[Tcur]+1,by=1,length.out=Tproj),
    rep('projected',Tproj),
    format(round(Cproj_region,1),nsmall=1),
    #c("C_2_5","C_25","C_50","C_75","C_97_5")
    c("C_025","C_25","C_50","C_75","C_975")
)
writemergedresults(df, 'Cproj_region', row.names=FALSE, quote=FALSE)
message('done Cproj_region')
# Cpred_region posterior predictive
Cpred_region_samples = load_samples('Cpred_region')
Cpred_region = t(apply(Cpred_region_samples,2,quantile,
    probs=c(.025,.25,.5,.75,.975)
))
df <- area_date_dataframe(
    quoted_regions,
    seq(dates[Tcond]+1,by=1,length.out=Mstep*Tstep),
    rep('inferred',Tlik),
    format(round(Cpred_region,1),nsmall=1),
    #c("C_2_5","C_25","C_50","C_75","C_97_5")
    c("C_025","C_25","C_50","C_75","C_975")
)
writemergedresults(df, 'Cpred_region', row.names=FALSE, quote=FALSE)
message('done Cpred_region')
})
}

##########################################################################
##########################################################################

epimap_cmdline_options = function(opt = Rmap_options()) {
  list(
    make_option(
      c("-s", "--spatialkernel"),
      type="character",
      default=opt$spatialkernel,
      help=paste(
          "Use spatial kernel (matern12/matern32/matern52/exp_quad/none);",
          "default =", opt$spatialkernel
      )
    ),
    make_option(
      c("-p", "--temporalkernel"),
      type="character",
      default=opt$temporalkernel,
      help=paste(
          "Use temporal kernel (matern12/matern32/matern52/exp_quad/none);",
          "default =", opt$temporalkernel
      )
    ),
    make_option(
      c("--gp_space_scale"),
      type="double",
      default=opt$gp_space_scale,
      help=paste(
          "If given and positive, set minimum space length scale to the value;",
          "default =", opt$gp_space_scale
      )
    ),

    make_option(
      c("--gp_time_scale"),
      type="double",
      default=opt$gp_time_scale,
      help=paste(
          "If given and positive, set minimum time length scale to the value;",
          "default =", opt$gp_time_scale
      )
    ),


    make_option(
      c("--fixed_gp_space_length_scale"),
      type="double",
      default=opt$fixed_gp_space_length_scale,
      help=paste(
          "If given and positive, fix the space length scale to the value;",
          "default =", opt$fixed_gp_space_length_scale
      )
    ),


    make_option(
      c("--fixed_gp_time_length_scale"),
      type="double",
      default=opt$fixed_gp_time_length_scale,
      help=paste(
          "If given and positive, fix the time length scale to the value;",
          "default =", opt$fixed_gp_time_length_scale
      )
    ),

    make_option(
      c("-l", "--localkernel"),
      type="character",
      default=opt$localkernel,
      help=paste("Use local kernel (local/none); default =", opt$localkernel)
    ),
    make_option(
      c("-g", "--globalkernel"),
      type="character",
      default=opt$globalkernel,
      help=paste("Use global kernel (global/none); default =", opt$globalkernel)
    ),
    make_option(
      c("-m", "--metapop"),
      type="character",
      default=opt$metapop,
      help=paste(
          "metapopulation model for inter-region cross infections",
          "(none, or comma separated list containing radiation{1,2,3},traffic{forward,reverse},uniform,in,in_out);",
          "default = ", opt$metapop
      )
    ),
    make_option(
      c("--constant_forward_rt"),  
      type="integer",
      default=opt$constant_forward_rt,
      help=paste("Use a the Rt from the last modelled week to predict forward; default =", opt$constant_forward_rt)
    ),
    make_option(
      c("--full_cases_distribution"),  
      type="integer",
      default=opt$constant_forward_rt,
      help=paste("Return the full distribution of cases, not just the distribution of the mean; default =", opt$full_cases_distribution)
    ),
    make_option(
      c("-v", "--observation_data"),
      type="character",
      default=opt$observation_data,
      help=paste(
          "observation values to use in the model",
          "(counts/clatent_mean/clatent_sample/clatent_recon/latent_reports);",
          "default =", opt$observation_data
      )
    ),
    make_option(
      c("-o", "--observation_model"),
      type="character",
      default=opt$observation_model,
      help=paste(
          "observation model",
          "(poisson/neg_binomial_{2,3}/gaussian);",
          "default =", opt$observation_model
      )
    ),
    make_option(
      c("-x", "--cleaned_sample_id"),
      type="integer",
      default=opt$cleaned_sample_id,
      help=paste("id of cleaned sample to use; default =", opt$cleaned_sample_id)
    ),
    make_option(
      c("-c", "--chains"),
      type="integer",
      default=opt$chains,
      help=paste("number of MCMC chains; default =", opt$chains)
    ),
    make_option(
      c("-i", "--iterations"),
      type="integer",
      default=opt$iterations,
      help=paste("Length of MCMC chains; defualt =", opt$iterations)
    ),
    make_option(
      c("--first_day_modelled"),
      type="character",
      default=opt$first_day_modelled,
      help=paste("Date of first day to model; default =",opt$first_day_modelled)
    ),
    make_option(
      c("--weeks_modelled"),
      type="integer",
      default=opt$weeks_modelled,
      help=paste("Number of weeks to model; default =",opt$weeks_modelled)
    ),
    make_option(
      c("--last_day_modelled"),
      type="character",
      default=opt$last_day_modelled,
      help=paste("Date of last day to model; default =",opt$last_day_modelled)
    ),
    make_option(
      c("--days_ignored"),
      type="integer",
      default=opt$days_ignored,
      help=paste("Days ignored; default =",opt$days_ignored)
    ),
    make_option(
      c("--days_predicted"),
      type="integer",
      default=opt$days_predicted,
      help=paste("Days predicted; default =",opt$days_predicted)
    ),
    make_option(
      c("--num_steps_forecasted"),
      type="integer",
      default=opt$num_steps_forecasted,
      help=paste("Days predicted; default =",opt$num_steps_forecasted)
    ),
    make_option(
      c("-d", "--results_directory"),
      type="character",
      default=opt$results_directory,
      help="If specified, store outputs in directory, otherwise use a unique directory"
    ),
    make_option(
      c("-r", "--clean_directory"),
      type="character",
      default=opt$clean_directory,
      help="Directory from which to load cleaned epidemic samples"
    ),
    make_option(
      c("--data_directory"),
      type="character",
      default=opt$data_directory,
      help="Directory from which to load data about regions and raw cases"
    ),
    make_option(
      c("-t", "--task_id"),
      type="integer",
      default=0,
      help="Task ID for Slurm usage. By default, turned off [0]."
    ),
    make_option(
      c("--limit_area"), 
      type="character", 
      default=opt$limit_area, 
      help=paste("If not NULL, center the radius of regions on this region",opt$limit_area)
    ),
    make_option(
      c("--limit_radius"), 
      type="double", 
      default=opt$limit_radius, 
      help=paste("If not NULL, the radius of regions to limit the data to; default",opt$limit_radius)
    )
  )
}

epimap_get_cmdline_options = function(opt=Rmap_options()) {
  cmdline_opt = epimap_cmdline_options(opt)
  opt_parser = OptionParser(option_list=cmdline_opt)
  parsed_opt = parse_args(opt_parser)
  for (name in names(parsed_opt)){
    opt[name] = parsed_opt[name]
  }
  opt
}
