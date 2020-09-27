library(rstan)
library(geosphere)
library(optparse)

Rmap_options = function(
  spatialkernel        = "matern12",
  temporalkernel       = "matern12",
  localkernel          = "local",
  globalkernel         = "global",
  metapop              = "radiation2,uniform,in",
  observation          = "cleaned_recon_sample",
  cleaned_sample_id    = 0,
  chains               = 1,
  iterations           = 6000,
  time_steps           = 15,
  days_per_step        = 7,
  days_ignored         = 6,
  days_predicted       = 2,
  num_steps_forecasted = 3,
  data_directory       = "data/",
  clean_directory      = "results/default",
  results_directory    = NULL
) {
     
  as.list(environment())
}

##########################################################################
##########################################################################
Rmap_read_data = function(env) { 
  with(env,{

    readdata = function(filename,...) {
      read.csv(sprintf('%s/%s.csv',opt$data_directory,filename),...)
    }
    readclean = function(filename,...) {
      read.csv(sprintf('%s/%s.csv',opt$clean_directory,filename),...)
    }

    #########################################################
    infprofile <- readdata("serial_interval")$fit
    Tip <- 30
    infprofile <- infprofile[1:Tip]
    infprofile <- infprofile/sum(infprofile)
    D <- length(infprofile)

    Tdp <- 14
    Adp <- 5.0
    Bdp <- 1.0
    testdelayprofile <- pgamma(1:Tdp,shape=Adp,rate=Bdp)
    testdelayprofile <- testdelayprofile/testdelayprofile[Tdp]
    testdelayprofile <- testdelayprofile - c(0.0,testdelayprofile[1:(Tdp-1)])


    df <- readdata("areas",row.names=1)
    N <- nrow(df)
    areas <- rownames(df)
    quoted_areas <- sapply(areas,function(s) paste('"',s,'"',sep=''))
    geoloc <- df[,1:2]
    population <- df[,3]

    geodist <- readdata("distances",row.names=1)
    colnames(geodist) <- areas

    # Use counts from uk_cases in case updated
    uk_cases <- readdata("uk_cases")
    ind <- sapply(uk_cases[,2], function(s)
        !(s %in% c('Outside Wales','Unknown','...17','...18'))
    )
    uk_cases <- uk_cases[ind,]
    AllCount <- uk_cases[,3:ncol(uk_cases)]
    Tall <- ncol(AllCount)
    dates <- as.Date(colnames(AllCount), format='X%Y.%m.%d')
    colnames(AllCount) <- dates
    rownames(AllCount) <- areas

    #########################################################
    radiation_length_scales <- c(.1,.2,.5)
    radiation_flux <- array(0,dim=c(N,N,length(radiation_length_scales)))
    for (i in 1:length(radiation_length_scales)) {
      ls <- radiation_length_scales[i]
      df <- data.matrix(readdata(sprintf('radiation_flux_ls=%1.1f',ls),row.names=1))
      radiation_flux[,,i] <- df
    }
    colnames(radiation_flux) <- areas
    rownames(radiation_flux) <- areas
    dimnames(radiation_flux)[[3]] <- radiation_length_scales

    traffic_flux <- array(0, dim=c(N,N,3))
    df <- data.matrix(readdata('traffic_flux_row-normed', row.names=1))
    traffic_flux[,,1] <- df
    df <- data.matrix(readdata('traffic_flux_max-normed', row.names=1))
    traffic_flux[,,2] <- df
    df <- data.matrix(readdata('traffic_flux_transpose_row-normed', row.names=1))
    traffic_flux[,,3] <- df
    colnames(traffic_flux) <- areas
    rownames(traffic_flux) <- areas

  })
  env
}# Rmap_read_data
##########################################################################
##########################################################################


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
    if (opt$observation == 'cleaned_latent_sample' ||
        opt$observation == 'cleaned_recon_sample') {
      sample_id = opt$cleaned_sample_id
      Clean_latent <- readclean(paste('Clatent_sample',sample_id,sep=''))
      Clean_recon <- readclean(paste('Crecon_sample',sample_id,sep=''))
      print(paste('Using samples from Clatent_sample',sample_id,'.csv',sep=''))
    } else if (opt$observation == 'cleaned_recon_sample') {
      sample_id = opt$cleaned_sample_id
      print(paste('Using samples from Crecon_sample',sample_id,'.csv',sep=''))
    } else {
      sample_id = 'mean'
      Clean_latent <- readclean('Clatent_mean')
      Clean_recon <- readclean('data/Crecon_median')
      # placeholder if not using cleaned data
    }

    #########################################################
    Mstep <- opt$time_steps        # Testing with 1 time period
    Tignore <- opt$days_ignored  # counts in most recent 7 days may not be reliable?
    if (!(Tall == length(Clean_latent) && Tall == length(Clean_recon))){
      print("WARNING: length of case data and cleaned data do not match. May need to regenerate the cleaned data. Truncating the case data")
    }
    Tall <- min(Tall, length(Clean_latent),length(Clean_recon))
    
    Tpred <- opt$days_predicted    # number of days held out for predictive probs eval
    Tstep <- opt$days_per_step # number of days to step for each time step of Rt prediction
    Tlik <- Mstep*Tstep     # number of days for likelihood to infer Rt
    Tall <- Tall-Tignore  # number of days; last 7 days counts ignore; not reliable
    Tcur <- Tall-Tpred       # number of days we condition on
    Tcond <- Tcur-Tlik       # number of days we condition on
    Mproj <- opt$num_steps_forecasted
    Tproj <- Tstep*Mproj           # number of days to project forward


    Count <- AllCount
    Count <- Count[,1:Tall] # get rid of ignored last days
    Clean_latent <- Clean_latent[,1:Tall] # get rid of ignored last days
    Clean_recon <- Clean_recon[,1:Tall] # get rid of ignored last days

    days_likelihood = seq(dates[Tcond+1],by=1,length.out=Tstep*Mstep)
    days_pred_held_out = seq(dates[Tcur+1],by=1,length.out=Tpred)

    print("Days used for likelihood fitting")
    print(days_likelihood)
    print("Days used for held out likelihood")
    print(days_pred_held_out)

    if (is.null(opt$results_directory)) {
      opt$results_directory = paste(
        'results/',
        as.character(Sys.time(),format='%Y%m%d'), 
        '-',as.character(days_likelihood[length(days_likelihood)],format='%Y%m%d'), 
        '-',opt$spatialkernel,  
        '-',opt$temporalkernel,  
        '-',opt$localkernel,  
        '-',opt$globalkernel,  
        '-',opt$metapop,  
        '-',opt$observation,
        sep=''
      )
    }
    dir.create(opt$results_directory, showWarnings = FALSE)
    if (opt$cleaned_sample_id>0) {
      runname = paste(opt$results_directory,'/',opt$cleaned_sample_id,'_',sep='')
    } else {
      runname = paste(opt$results_directory,'/',sep='')
    }
    print(runname)

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
    OBSERVATIONMODELS = list(
      'poisson' = 1,
      'neg_binomial_2' = 2,
      'neg_binomial_3' = 3,
      'cleaned_latent_mean' = 4,
      'cleaned_latent_sample' = 4,
      'cleaned_recon_sample' = 5
    )
    OBSERVATIONMODEL = OBSERVATIONMODELS[[opt$observation]]
    if (is.null(OBSERVATIONMODEL)) {
      stop(c('Unrecognised observation option ',opt$observation));
    }

    #########################################################
    # metapopulation cross-area fluxes.
    METAPOPMODEL = strsplit(opt$metapop,',')[[1]]
    METAPOPOPTIONS = list(
      'radiation1' = radiation_flux[,,1], # smoothed radiation model with length scale = .1 (10km)
      'radiation2' = radiation_flux[,,2], # smoothed radiation model with length scale = .2 (20km)
      'radiation3' = radiation_flux[,,3], # smoothed radiation model with length scale = .5 (50km)
      'traffic1' = traffic_flux[,,1], # traffic flux (all row sums = 1)
      'traffic2' = traffic_flux[,,2], # traffic flux (max row sum = 1)
      'traffic3' = traffic_flux[,,3], # traffic flux transpose (all row sums = 1) 
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
    times = 1:Mstep
    timedist = matrix(0, Mstep, Mstep)
    for (i in 1:Mstep) {
      for (j in 1:Mstep) {
        timedist[i, j] = abs(times[i] - times[j]) * Tstep
      }
    }

    # precompute lockdown cutoff kernel
    lockdown_day = as.Date("2020-03-23")
    days_lik_start = days_likelihood[seq(1, length(days_likelihood), Tstep)]
    days_lik_start = vapply(days_lik_start, (function (day) as.Date(day, format="%Y-%m-%d")), double(1))
    day_pre_lockdown = vapply(days_lik_start, (function (day) day < lockdown_day), logical(1))
    
    time_corellation_cutoff = matrix(0,Mstep,Mstep)
    for (i in 1:Mstep) {
      for (j in 1:Mstep) {
        time_corellation_cutoff[i, j] = !xor(day_pre_lockdown[i], day_pre_lockdown[j])
      }
    }

    #########################################################
    # Main computation
    Rmap_data <- list(
      N = N, 
      Mstep = Mstep,
      Tall = Tall,
      Tcond = Tcond,
      Tlik = Tlik,
      Tproj = Tproj,
      Tstep=Tstep,

      Count = Count,
      Clean_latent = Clean_latent,
      Clean_recon = Clean_recon,
      geodist = geodist,
      timedist = timedist,
      timecorcut = time_corellation_cutoff,

      SPATIAL_KERNEL = SPATIAL_KERNEL,
      TEMPORAL_KERNEL = TEMPORAL_KERNEL,
      LOCAL_KERNEL = LOCAL_KERNEL,
      GLOBAL_KERNEL = GLOBAL_KERNEL,
      DO_METAPOP = DO_METAPOP,
      DO_IN_OUT = DO_IN_OUT,
      OBSERVATIONMODEL = OBSERVATIONMODEL,

      Tip = Tip, 
      infprofile = infprofile,
      Tdp = Tdp,
      delayprofile = testdelayprofile,
      F = F,
      flux = flux
    )

    #########################################################
    Rmap_init = list()
    for (i in 1:numchains) {
      Rmap_init[[i]] = list(
        gp_time_length_scale = 100.0,
        gp_space_length_scale = 2.0,
        gp_space_sigma = .01,
        global_sigma = .01,
        local_scale = .01,
        dispersion = 5.0
      )
    }


    #########################################################
    fit <- stan(file = 'stan_files/Rmap.stan',
                data = Rmap_data, 
                init = Rmap_init,
                iter = numiters, 
                chains = numchains,
                control = list(adapt_delta = .9))


    #########################################################
    saveRDS(fit, paste(runname,'stanfit.rds', sep=''))

    #########################################################
    # Summary of fit
    print(summary(fit, 
        pars=c("gp_space_length_scale","gp_space_sigma","gp_time_length_scale",
            "global_sigma","local_scale","dispersion",
            "Rt_all","coupling_rate","flux_probs"), 
        probs=0.5)$summary)


    #########################################################
    end_time <- Sys.time()
    print("Time to run")
    print(end_time - start_time)
    print(runname)

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
    save_samples = function(pars,areafirst=FALSE) {
      samples = extract(fit,pars=pars,permuted=FALSE)
      samples = samples[seq(numchains,by=numchains,to=numiters/2),,,drop=FALSE]
      ns = numiters/2
      np = dim(samples)[3]/N
      parnames = dimnames(samples)[[3]]
      dim(samples) = c(ns,np*N)
      colnames(samples) = parnames
      if (areafirst) {
        ind = N*(rep(1:np,N)-1) + rep(1:N,each=np) 
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
    
    times = rep(c(1:Mstep,rep(Mstep,Mproj)), N)
    places = rep(1:N, each=Mstep+Mproj)
    indicies = places + (N)*(times-1)
    Rt = Rt[indicies,]

    Rt <- Rt[sapply(1:(N*(Mstep+Mproj)),function(i)rep(i,Tstep)),]
    print(sprintf("median Rt range: [%f, %f]",min(Rt[,"50%"]),max(Rt[,"50%"])))
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
    dim(Rt) <- c(numsamples,Mstep,N)
    Pexceedance = array(0.0,dim=c(Mstep,N,numthresholds))
    for (k in 1:Mstep) {
      for (i in 1:N) {
        for (x in 1:numthresholds) {
          Pexceedance[k,i,x] = mean(Rt[,k,i]>thresholds[x])
        }
      }
    }
    Pexceedance = Pexceedance[c(1:Mstep,rep(Mstep,Mproj)),,]
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
    print(sprintf("median Cpred range: [%f, %f]",min(Cpred[,"50%"]),max(Cpred[,"50%"])))
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
    Cweekly <- as.matrix(Count[,(Tcond+1):(Tcond+Tlik)])
    dim(Cweekly) <- c(N,Tstep,Mstep)
    Cweekly <- apply(Cweekly,c(1,3),sum)

    ignoredweek <- apply(AllCount[,(Tcur+1):(Tcur+Tpred+Tignore)],c(1),sum)
    Cweekly <- cbind(Cweekly,ignoredweek)

    s <- summary(fit, pars=c("Cproj"), probs=c(.5))$summary
    projectedweeks <- as.matrix(s[,"50%"])
    dim(projectedweeks) <- c(Tstep,Mproj,N)
    projectedweeks <- projectedweeks[,2:Mproj,,drop=FALSE]
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
      format(round(Cweekly,1),digits=6),
      c("C_weekly")
    )
    writeresults(df, 'Cweekly', row.names=FALSE, quote=FALSE)



    s <- summary(fit, pars=c("Cproj"), probs=c(0.025, .25, .5, .75, .975))$summary
    Cproj <- s[,c("2.5%","25%", "50%","75%", "97.5%")]
    #Cproj <- s[,c("2.5%", "50%", "97.5%")]
    Cproj <- t(t(Cproj))
    print(sprintf("median Cproj range: [%f, %f]",min(Cproj[,"50%"]),max(Cproj[,"50%"])))
    df <- area_date_dataframe(
      quoted_areas,
      seq(dates[Tcur]+1,by=1,length.out=Tproj),
      rep('projected',Tproj),
      format(round(Cproj,1),digits=5),
      #c("C_2_5","C_25","C_50","C_75","C_97_5")
      c("C_025","C_25","C_50","C_75","C_975")
    )
    writeresults(df, 'Cproj.csv', row.names=FALSE, quote=FALSE)


    #################################################################
    # predictive probabilities
    s <- summary(fit, pars="Ppred", probs=c(0.025, .5, .975))$summary
    Ppred <- s[,"mean"]
    logpred <- log(Ppred)
    dim(logpred) <- c(Tpred,N)
    logpred <- t(logpred)
    print(sprintf("mean log predictives = %f",mean(logpred)))
    df <- data.frame(area = quoted_areas, logpred = logpred, provenance=rep('inferred', length(quoted_areas)))
    for (i in 1:Tpred)
      colnames(df)[i+1] <- sprintf('logpred_day%d',i)
    writeresults(df, 'logpred', row.names=FALSE, quote=FALSE)


    ####################################################################
    # pairs plot
    #pdf(paste('fits/',runname,'_pairs.pdf',sep=''),width=9,height=9)
    #pairs(fit, pars=c(
    #    "gp_space_length_scale","gp_space_sigma","gp_time_length_scale",
    #    "global_sigma","local_scale","precision")) 
    #dev.off()

    ####################################################################
  })
  env
} # Rmap_postprocess
##########################################################################
##########################################################################


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

    Rt_samples = load_samples('Rt')
    Cpred_samples = load_samples('Cpred')
    Cproj_samples = load_samples('Cproj')

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
    Rt = t(apply(Rt_samples,2,quantile,
        probs=c(0.025, .1, .2, 0.25, .3, .4, .5, .6, .7, 0.75, .8, .9, .975)
    ))
    
    Rt = Rt[sapply(1:N,function(i)rep((i-1)*Mstep+c(1:Mstep,rep(Mstep,Mproj)),each=Tstep)),]
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

    #################################################################
    # Rt exceedance probabilities
    thresholds = c(.8, .9, 1.0, 1.1, 1.2, 1.5, 2.0)
    numthresholds = length(thresholds)
    numsamples = numruns * numiters/2
    Rt <- as.matrix(Rt_samples)
    dim(Rt) <- c(numsamples,Mstep,N)
    Pexceedance = array(0.0,dim=c(Mstep,N,numthresholds))
    for (k in 1:Mstep) {
      for (i in 1:N) {
        for (x in 1:numthresholds) {
          Pexceedance[k,i,x] = mean(Rt[,k,i]>thresholds[x])
        }
      }
    }
    Pexceedance = Pexceedance[c(1:Mstep,rep(Mstep,Mproj)),,]
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



    #################################################################
    # posterior predictives and projections
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

    ####################################################################################
    # weekly counts. Includes 1 last column of actual counts among days ignored in model
    Cweekly <- as.matrix(Count[,(Tcond+1):(Tcond+Tlik)])
    dim(Cweekly) <- c(N,Tstep,Mstep)
    Cweekly <- apply(Cweekly,c(1,3),sum)
    
    ignoredweek <- apply(AllCount[,(Tcur+1):(Tcur+Tpred+Tignore)],c(1),sum)
    Cweekly <- cbind(Cweekly,ignoredweek)

    projectedweeks = as.matrix(apply(Cproj_samples,2,quantile,
        probs=c(.5)
    ))
    dim(projectedweeks) <- c(Tstep,Mproj,N)
    projectedweeks <- projectedweeks[,2:Mproj,,drop=FALSE]
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



    Cproj = t(apply(Cproj_samples,2,quantile,
        probs=c(.025,.25,.5,.75,.975)
    ))
    df <- area_date_dataframe(
        quoted_areas,
        seq(dates[Tcur]+1,by=1,length.out=Tproj),
        rep('projected',Tproj),
        format(round(Cproj,1),digits=1),
        #c("C_2_5","C_25","C_50","C_75","C_97_5")
        c("C_025","C_25","C_50","C_75","C_975")
    )
    writemergedresults(df, 'Cproj', row.names=FALSE, quote=FALSE)

  })
}
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

Rt_samples = load_samples('Rt')
Cpred_samples = load_samples('Cpred')
Cproj_samples = load_samples('Cproj')

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
Rt = t(apply(Rt_samples,2,quantile,
    probs=c(0.025, .1, .2, 0.25, .3, .4, .5, .6, .7, 0.75, .8, .9, .975)
))

Rt = Rt[sapply(1:N,function(i)rep((i-1)*Mstep+c(1:Mstep,rep(Mstep,Mproj)),each=Tstep)),]
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

#################################################################
# Rt exceedance probabilities
thresholds = c(.8, .9, 1.0, 1.1, 1.2, 1.5, 2.0)
numthresholds = length(thresholds)
numsamples = numruns * numiters/2
Rt <- as.matrix(Rt_samples)
dim(Rt) <- c(numsamples,Mstep,N)
Pexceedance = array(0.0,dim=c(Mstep,N,numthresholds))
for (k in 1:Mstep) {
  for (i in 1:N) {
    for (x in 1:numthresholds) {
      Pexceedance[k,i,x] = mean(Rt[,k,i]>thresholds[x])
    }
  }
}
Pexceedance = Pexceedance[c(1:Mstep,rep(Mstep,Mproj)),,]
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



#################################################################
# posterior predictives and projections
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

####################################################################################
# weekly counts. Includes 1 last column of actual counts among days ignored in model
Cweekly <- as.matrix(Count[,(Tcond+1):(Tcond+Tlik)])
dim(Cweekly) <- c(N,Tstep,Mstep)
Cweekly <- apply(Cweekly,c(1,3),sum)

ignoredweek <- apply(AllCount[,(Tcur+1):(Tcur+Tpred+Tignore)],c(1),sum)
Cweekly <- cbind(Cweekly,ignoredweek)

projectedweeks = as.matrix(apply(Cproj_samples,2,quantile,
    probs=c(.5)
))
dim(projectedweeks) <- c(Tstep,Mproj,N)
projectedweeks <- projectedweeks[,2:Mproj,,drop=FALSE]
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



Cproj = t(apply(Cproj_samples,2,quantile,
    probs=c(.025,.25,.5,.75,.975)
))
df <- area_date_dataframe(
    quoted_areas,
    seq(dates[Tcur]+1,by=1,length.out=Tproj),
    rep('projected',Tproj),
    format(round(Cproj,1),digits=1),
    #c("C_2_5","C_25","C_50","C_75","C_97_5")
    c("C_025","C_25","C_50","C_75","C_975")
)
writemergedresults(df, 'Cproj', row.names=FALSE, quote=FALSE)

})
}

##########################################################################
##########################################################################
