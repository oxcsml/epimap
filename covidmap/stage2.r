library(rstan)
library(geosphere)
library(optparse)
source("covidmap/read_data.r")
source('covidmap/utils.r')
source("epimap/epimap.r")

covidmap_stage2_options = function(
  spatialkernel        = "matern12",
  temporalkernel       = "matern12",
  localkernel          = "local",
  globalkernel         = "global",
  gp_space_scale       = 0.5, # units of 100km
  gp_space_decay_scale = .25,
  gp_time_scale        = 50.0, # units of 1 day
  gp_time_decay_scale  = .25,
  fixed_gp_space_length_scale = .1, # 0.5,
  fixed_gp_time_length_scale = 100, # 100.0,
  constant_forward_rt  = 0, # if 0, use the Rt from last week to predict forwards
  full_cases_distribution = 1,
  metapop              = "in,alt_traffic_forward,alt_traffic_reverse",
  #metapop              = "traffic_forward,traffic_reverse,radiation1,radiation2,radiation3,uniform,in",
  approximation        = "twostage",
  observation_data     = "cleaned_latent_sample",
  observation_model    = "gaussian",
  singlearea_sample_id    = 0, 
  region_id            = 0,

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
  singlearea_directory = "fits/singlearea", 
  results_directory    = NULL,

  limit_area           = NULL,
  limit_radius         = NULL

) {
  as.list(environment())
}

##########################################################################
##########################################################################

covidmap_stage2_setup = function(opt = covidmap_stage2_options()) {
  env = new.env(parent=globalenv())
  env$opt = opt
  covidmap_read_data(env)
  with(env,{

    stopifnot(identical(opt$approximation,"twostage") || identical(opt$approximation,"regional"))

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
    stopifnot(Tcur+Tproj == length(Clean_latent))
    message("Tcur: ", Tcur, ", length Clean_recon: ", length(Clean_recon))
    stopifnot(Tcur == length(Clean_recon))
    print(length(Clean_latent))

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
    if (opt$singlearea_sample_id>0 || opt$region_id>0) {
      runname = paste(
        opt$results_directory,'/',
        opt$region_id,'_',
        opt$singlearea_sample_id,'_',
        sep=''
      )
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

covidmap_stage2_run = function(
    env
) {
    with(env, {
        start_time <- Sys.time()
    
        time_dists = time_distances(
            Mstep + Mproj,
            Tstep,
            days,
        )

        metapop_list = strsplit(opt$metapop,',')[[1]]
        metapop_type = metapop_list[1]
        metapop_list = metapop_list[2:length(metapop_list)]
        flux_options = list(
            # 'radiation1' = radiation_flux[,,1], # smoothed radiation model with length scale = .1 (10km)
            # 'radiation2' = radiation_flux[,,2], # smoothed radiation model with length scale = .2 (20km)
            # 'radiation3' = radiation_flux[,,3], # smoothed radiation model with length scale = .5 (50km)
            'traffic_forward' = traffic_flux[,,1], # infected commuters taking infection from home to work
            'traffic_reverse' = traffic_flux[,,2], # commuters getting infected at work and bringing back home
            'alt_traffic_forward' = alt_traffic_flux[,,1], # infected commuters taking infection from home to work
            'alt_traffic_reverse' = alt_traffic_flux[,,2], # commuters getting infected at work and bringing back home
            'uniform' = matrix(1.0/N,N,N) # uniform cross-area infection
        )
        if (length(metapop_list) > 1) {
            flux = list()
            F = 0
            for (s in metapop_list) {
                t = flux_options[[s]]
                if (is.null(t)) {
                    stop(c('Unrecognised metapop option ',s));
                }
                F = F+1
                flux[[F]] = t
            }
        } else {
            flux = NULL
        }

        if (identical(opt$approximation, "regional")) {
            fit <- epimap_region(
                iter = numiters,
                chains = numchains,
                cases = AllCount,
                latent_infections = Clean_latent,
                Tcur = Tcur,
                Tstep = Tstep,
                Tpred = Tpred,
                Mstep = Mstep,
                Mignore = opt$steps_ignored_stage2,
                Mproj = Mproj,
                area_modelled = modelled_region[,opt$region_id],
                area_inferred = inferred_region[,opt$region_id],
                spatialkernel = opt$spatialkernel,
                temporalkernel = opt$temporalkernel,
                localkernel = opt$localkernel,
                globalkernel = opt$globalkernel,
                observation_data = opt$observation_data,
                observation_model = opt$observation_model,
                infection_profile = infprofile,
                delay_profile = testdelayprofile,
                geographic_distances = geodist,
                temporal_distances = time_dists$time_distances,
                temporal_correlation = time_dists$time_correlations,
                metapop = metapop_type,
                flux = flux,
                gp_space_scale = opt$gp_space_scale,
                gp_space_decay_scale = opt$gp_space_decay_scale,
                gp_time_scale = opt$gp_time_scale,
                gp_time_decay_scale = opt$gp_time_decay_scale,
                fixed_gp_space_length_scale = opt$fixed_gp_space_length_scale,
                fixed_gp_time_length_scale = opt$fixed_gp_time_length_scale,
                constant_forward_rt = opt$constant_forward_rt,
                full_cases_dist = opt$full_cases_distribution
            )
        } else if (identical(opt$approximation, "twostage")) {
            fit <- epimap_twostage(
                iter = numiters,
                chains = numchains,
                cases = AllCount,
                reconstructed_infections = Clean_recon,
                latent_infections = Clean_latent,
                Tcur = Tcur,
                Tstep = Tstep,
                Tpred = Tpred,
                Mstep = Mstep,
                Mignore = opt$steps_ignored_stage2,
                Mproj = Mproj,
                spatialkernel = opt$spatialkernel,
                temporalkernel = opt$temporalkernel,
                localkernel = opt$localkernel,
                globalkernel = opt$globalkernel,
                observation_data = opt$observation_data,
                observation_model = opt$observation_model,
                sparse_region = sparse_region,
                infection_profile = infprofile,
                delay_profile = testdelayprofile,
                geographic_distances = geodist,
                temporal_distances = time_dists$time_distances,
                temporal_correlation = time_dists$time_correlations,
                metapop = metapop_type,
                flux = flux,
                gp_space_scale = opt$gp_space_scale,
                gp_space_decay_scale = opt$gp_space_decay_scale,
                gp_time_scale = opt$gp_time_scale,
                gp_time_decay_scale = opt$gp_time_decay_scale,
                fixed_gp_space_length_scale = opt$fixed_gp_space_length_scale,
                fixed_gp_time_length_scale = opt$fixed_gp_time_length_scale,
                constant_forward_rt = opt$constant_forward_rt,
                full_cases_dist = opt$full_cases_distribution
            )
        } else {
            stop(cat("Unrecognised approximation: ", opt$approximation))
        }

        saveRDS(fit, paste(runname,'stanfit.rds', sep=''))

        end_time <- Sys.time()
        message(paste(
        "Time to run:",
        round(difftime(end_time, start_time, units="hours"),3),
        "hours")
        )
        message(runname)
    })
    env
}

##########################################################################
##########################################################################

covidmap_stage2_postprocess = function(env) {
  with(env, {
    writeresults = function(data,filename,...) {
      write.csv(data,sprintf('%s%s.csv',runname,filename),...)
    }

    if (opt$region_id>0) {
      Ninferred = sum(inferred_region[,opt$region_id])
      inferred_areas = inferred_region[,opt$region_id]
    } else {
      Ninferred = N
      inferred_areas = 1:N
    }

    # save raw samples
    save_samples = function(pars,N_sites=Ninferred,areafirst=FALSE) {
      samples = extract(fit,pars=pars,permuted=FALSE)
      samples = samples[seq(opt$thinning,by=opt$thinning,to=numiters/2),,,drop=FALSE]
      ns = numiters/2/opt$thinning
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
    if (opt$region_id>0) {
      nregion = 1
    } else {
      nregion = N_region
    }
    save_samples("Rt_region", N_sites=nregion)
    save_samples("Cpred_region", N_sites=nregion, areafirst=TRUE)
    save_samples("Cproj_region", N_sites=nregion, areafirst=TRUE)

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
    print("Rt")
    s <- summary(fit, pars="Rt", probs=c(0.025, .1, .2, 0.25, .3, .4, .5, .6, .7, 0.75, .8, .9, .975))$summary
    Rt <- s[,c("2.5%","10%", "20%", "25%", "30%", "40%", "50%", "60%", "70%", "75%", "80%", "90%","97.5%")]
    #s <- summary(fit, pars="Rt", probs=c(.1, .2, .3, .4, .5, .6, .7, .8, .9))$summary
    #Rt <- s[,c("10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%")]
    
    times = rep(c(1:(Mstep+Mproj)), Ninferred)
    places = rep(1:Ninferred, each=Mstep+Mproj)
    indicies = places + (Ninferred)*(times-1)
    Rt = Rt[indicies,]

    Rt <- Rt[sapply(1:(Ninferred*(Mstep+Mproj)),function(i)rep(i,Tstep)),]
    message(sprintf("median Rt range: [%f, %f]",min(Rt[,"50%"]),max(Rt[,"50%"])))
    df <- area_date_dataframe(
      quoted_areas[inferred_areas],
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
    print("Pexceed_samples")
    thresholds = c(.8, .9, 1.0, 1.1, 1.2, 1.5, 2.0)
    numthresholds = length(thresholds)
    numsamples = numchains*numiters/2
    Rt <- as.matrix(fit, pars="Rt")
    dim(Rt) <- c(numsamples,Mstep+Mproj,Ninferred)
    Pexceedance = array(0.0,dim=c(Mstep+Mproj,Ninferred,numthresholds))
    for (k in 1:(Mstep+Mproj)) {
      for (i in 1:Ninferred) {
        for (x in 1:numthresholds) {
          Pexceedance[k,i,x] = mean(Rt[,k,i]>thresholds[x])
        }
      }
    }
    Pexceedance = Pexceedance[c(1:(Mstep+Mproj)),,]
    Pexceedance <- Pexceedance[sapply(1:(Mstep+Mproj),function(k)rep(k,Tstep)),,]
    dim(Pexceedance) <- c(Tstep*(Mstep+Mproj)*Ninferred,numthresholds)
    df <- area_date_dataframe(
        quoted_areas[inferred_areas],
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
      quoted_areas[inferred_areas],
      seq(dates[Tcond]+1,by=1,length.out=Mstep*Tstep),
      rep('inferred',Tlik),
      format(round(Cpred,1),nsmall=1),
      #c("C_2_5","C_25","C_50","C_75","C_97_5")
      c("C_025","C_25","C_50","C_75","C_975")
    )
    writeresults(df, 'Cpred', row.names=FALSE, quote=FALSE)


    # weekly counts. Includes 1 last column of actual counts among days ignored in model
    Tweek = Tstep # assumes Tstep = 7
    Cweekly <- as.matrix(AllCount[inferred_areas,(Tcond+1):(Tcond+Tlik)]) 
    dim(Cweekly) <- c(Ninferred,Tstep,Mstep)
    Cweekly <- apply(Cweekly,c(1,3),sum)

    # ignoredweek <- apply(AllCount[,(Tcur+1):(Tcur+Tstep)],c(1),sum)
    # Cweekly <- cbind(Cweekly,ignoredweek)

    s <- summary(fit, pars=c("Cproj"), probs=c(.5))$summary
    projectedweeks <- as.matrix(s[,"50%"])
    dim(projectedweeks) <- c(Tstep,Mproj,Ninferred)
    projectedweeks <- projectedweeks[,1:Mproj,,drop=FALSE]
    projectedweeks <- apply(projectedweeks,c(2,3),sum)
    projectedweeks <- t(projectedweeks)

    Cweekly <- cbind(Cweekly,projectedweeks)
    Cweekly <- t(Cweekly)
    Cweekly <- Cweekly[sapply(1:(Mstep+Mproj),function(k)rep(k,Tstep)),]
    dim(Cweekly) <- c(Ninferred*(Tlik+Tstep*Mproj))

    Cweekly_provenance <- c(rep('actual',Tlik),rep('projected',Tproj))
    df <- area_date_dataframe(
      quoted_areas[inferred_areas],
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
      quoted_areas[inferred_areas], 
      seq(dates[Tcur]+1,by=1,length.out=Tproj),
      rep('projected',Tproj),
      format(round(Cproj,1),digits=5),
      #c("C_2_5","C_25","C_50","C_75","C_97_5")
      c("C_025","C_25","C_50","C_75","C_975")
    )
    writeresults(df, 'Cproj', row.names=FALSE, quote=FALSE)


    #################################################################
    # predictive probabilities
    print("logpred")
    s <- summary(fit, pars="Ppred", probs=c(0.025, .5, .975))$summary
    Ppred <- s[,"mean"]
    logpred <- log(Ppred)
    dim(logpred) <- c(Tpred,Ninferred)
    logpred <- t(logpred)
    dim(logpred) = c(Ninferred*Tpred)
    message(sprintf("mean log predictives = %f",mean(logpred)))
    df <- area_date_dataframe(
      quoted_areas[inferred_areas],
      seq(dates[Tcur]+1,by=1,length.out=Tpred),
      rep('projected', Tpred),
      logpred,
      c('logpred')
    )
    writeresults(df, 'logpred', row.names=FALSE, quote=FALSE)
  })
  env
}

##########################################################################
##########################################################################

covidmap_stage2_merge = function(
    env,
    singlearea_sample_ids=c(0),
    region_ids=c(0)
) {
  env$singlearea_sample_ids = singlearea_sample_ids
  env$region_ids = region_ids
  with(env, {

  writemergedresults = function(data,filename,...) {
    write.csv(data,sprintf('%s/merged_%s.csv',opt$results_directory,filename),...)
  }

  numruns = length(singlearea_sample_ids)

  load_samples = function(region_id,pars) {
    samples = do.call(rbind, lapply(1:numruns, function(i) {
      read.csv(paste(
        opt$results_directory,
        '/',region_id,
        '_',singlearea_sample_ids[i],
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
  areas = function(region_id) {
    if (region_id==0) {
      rep(1,N)==1
    } else {
      inferred_region[,region_id]==1
    }
  }
  Nareas = function(region_id) sum(areas(region_id))

  if (length(region_ids)==1 && region_ids[1]==0) {
    quoted_areas_by_regions = quoted_areas
  } else {
    quoted_areas_by_regions = do.call(c,lapply(region_ids, function(region_id) {
      quoted_areas[areas(region_id)]
    }))
  }

  #################################################################
  # Rt posterior
  Rt = do.call(rbind,lapply(region_ids, function(region_id) {

    Rt_samples = load_samples(region_id,'Rt')
    # TODO: we get very infrequent Nas in the cori model
    # if (any(is.na(Rt_samples))) {
    #   message("WARNING: NAs in Rt samples")
    #   Rt_samples = Rt_samples[complete.cases(Rt_samples),]
    # }
    
    t(apply(Rt_samples,2,quantile,
      probs=c(0.025, .1, .2, 0.25, .3, .4, .5, .6, .7, 0.75, .8, .9, .975)
    ))
  }))

  Rt = Rt[sapply(1:N,function(i)rep((i-1)*(Mstep+Mproj)+c(1:(Mstep+Mproj)),each=Tstep)),]
  #Rt = Rt[sapply(1:N,function(i)rep((i-1)*Mstep+c(1:Mstep,rep(Mstep,Mproj)),each=Tstep)),]
  df <- area_date_dataframe(
      quoted_areas_by_regions,
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
  #numsamples = numruns * floor(numiters/2 / opt$thinning)

  Pexceedance = do.call(rbind,lapply(region_ids, function(region_id) {
    Rt_samples = load_samples(region_id,'Rt')
    Rt_samples <- as.matrix(Rt_samples)
    numsamples = dim(Rt_samples)[1] 
    Narea = Nareas(region_id)
    dim(Rt_samples) <- c(numsamples,Mstep+Mproj,Narea)
    Pexceedance = array(0.0,dim=c(Mstep+Mproj,Narea,numthresholds))
    for (k in 1:(Mstep+Mproj)) {
      for (i in 1:Narea) {
        for (x in 1:numthresholds) {
          Pexceedance[k,i,x] = mean(Rt_samples[,k,i]>thresholds[x])
        }
      }
    }
    Pexceedance = Pexceedance[c(1:(Mstep+Mproj)),,]
    Pexceedance <- Pexceedance[sapply(1:(Mstep+Mproj),function(k)rep(k,Tstep)),,]
    dim(Pexceedance) <- c(Tstep*(Mstep+Mproj)*Narea,numthresholds)
    Pexceedance
  }))
  df <- area_date_dataframe(
      quoted_areas_by_regions,
      days_all,
      provenance,
      format(round(Pexceedance,2),nsmall=2),
      c("P_08","P_09","P_10","P_11","P_12","P_15","P_20")
  )
  writemergedresults(df, 'Pexceed', row.names=FALSE, quote=FALSE)
  message('done Pexceedance')


  #################################################################
  # posterior predictives and projections
  Cpred = do.call(rbind,lapply(region_ids, function(region_id) {
    Cpred_samples = load_samples(region_id,'Cpred')
    t(apply(Cpred_samples,2,quantile,
      probs=c(0.025, 0.25, .5, 0.75, .975)
    ))
  }))

  df <- area_date_dataframe(
      quoted_areas_by_regions,
      seq(dates[Tcond]+1,by=1,length.out=Mstep*Tstep),
      rep('inferred',Tlik),
      format(round(Cpred,1),nsmall=1),
      #c("C_2_5","C_25","C_50","C_75","C_97_5")
      c("C_025","C_25","C_50","C_75","C_975")
  )
  writemergedresults(df, 'Cpred', row.names=FALSE, quote=FALSE)
  message('done Cpred')

  ####################################################################################
  # weekly counts. Includes 1 last column of actual counts among days ignored in model

  Cweekly <- as.matrix(AllCount[,(Tcond+1):(Tcond+Tlik)])
  dim(Cweekly) <- c(N,Tstep,Mstep)
  Cweekly <- apply(Cweekly,c(1,3),sum)

  stopifnot(Tcur+Tstep<=length(AllCount))

  #ignoredweek <- apply(AllCount[,(Tcur+1):(Tcur+Tstep)],c(1),sum)
  #Cweekly <- cbind(Cweekly,ignoredweek)

  Cweekly = do.call(rbind,lapply(region_ids,function(region_id) {
    Narea = Nareas(region_id)
    Cproj_samples = load_samples(region_id,'Cproj')
    projectedweeks = as.matrix(apply(Cproj_samples,2,quantile,
        probs=c(.5)
    ))
    dim(projectedweeks) <- c(Tstep,Mproj,Narea)
    projectedweeks <- projectedweeks[,1:Mproj,,drop=FALSE]
    projectedweeks <- apply(projectedweeks,c(2,3),sum)
    projectedweeks <- t(projectedweeks)
    cbind(Cweekly[areas(region_id),],projectedweeks)
  }))

  Cweekly <- t(Cweekly)
  Cweekly <- Cweekly[sapply(1:(Mstep+Mproj),function(k)rep(k,Tstep)),]
  dim(Cweekly) <- c(N*(Tlik+Tstep*Mproj))
  df <- area_date_dataframe(
      quoted_areas_by_regions,
      days_all,
      provenance,
      format(Cweekly,digits=3),
      c("C_weekly")
  )
  writemergedresults(df, 'Cweekly', row.names=FALSE, quote=FALSE)
  message('done Cweekly')



  Cproj = do.call(rbind,lapply(region_ids,function(region_id) {
    Cproj_samples = load_samples(region_id,'Cproj')
    t(apply(Cproj_samples,2,quantile,
        probs=c(.025,.25,.5,.75,.975)
    ))
  }))
  df <- area_date_dataframe(
      quoted_areas_by_regions,
      seq(dates[Tcur]+1,by=1,length.out=Tproj),
      rep('projected',Tproj),
      format(round(Cproj,1),nsmall=1),
      #c("C_2_5","C_25","C_50","C_75","C_97_5")
      c("C_025","C_25","C_50","C_75","C_975")
  )
  writemergedresults(df, 'Cproj', row.names=FALSE, quote=FALSE)
  message('done Cproj')


  # Rt_region posterior
  Rt_region = do.call(rbind,lapply(region_ids,function(region_id) {
    Rt_region_samples = load_samples(region_id,'Rt_region')
    t(apply(Rt_region_samples,2,quantile,
      probs=c(0.025, .1, .2, 0.25, .3, .4, .5, .6, .7, 0.75, .8, .9, .975)
    ))
  }))

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
  Cproj_region = do.call(rbind,lapply(region_ids,function(region_id) {
    Cproj_region_samples = load_samples(region_id,'Cproj_region')
    t(apply(Cproj_region_samples,2,quantile,
      probs=c(.025,.25,.5,.75,.975)
    ))
  }))
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
  Cpred_region = do.call(rbind,lapply(region_ids,function(region_id) {
    Cpred_region_samples = load_samples(region_id,'Cpred_region')
    t(apply(Cpred_region_samples,2,quantile,
      probs=c(.025,.25,.5,.75,.975)
    ))
  }))
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

covidmap_stage2_cmdline_options = function(opt = covidmap_stage2_options()) {
  list(
    make_option(
      c("--approximation"),
      type="character",
      default=opt$localkernel,
      help=paste("Run approximation (singlearea/twostage/regional); default =", opt$approximation)
    ),

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
          "(none, or comma separated list containing radiation{1,2,3},{alt_}traffic_{forward,reverse},uniform,in,in_out);",
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
      c("-x", "--singlearea_sample_id"),
      type="integer",
      default=opt$singlearea_sample_id,
      help=paste("id of cleaned sample to use; default =", opt$singlearea_sample_id)
    ),
    make_option(
      c("--region_id"),
      type="integer",
      default=opt$region_id,
      help=paste("id of region to model; default =", opt$region_id)
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
      c("-r", "--singlearea_directory"),
      type="character",
      default=opt$singlearea_directory,
      help="Directory from which to load single area approximation epidemic samples"
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

covidmap_stage2_get_cmdline_options = function(opt=covidmap_stage2_options()) {
  cmdline_opt = covidmap_stage2_cmdline_options(opt)
  opt_parser = OptionParser(option_list=cmdline_opt)
  parsed_opt = parse_args(opt_parser)
  for (name in names(parsed_opt)){
    opt[name] = parsed_opt[name]
  }
  opt
}
