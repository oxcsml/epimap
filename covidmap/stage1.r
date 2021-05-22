library(optparse)
library(gsubfn)
library(plyr)
source("covidmap/read_data.r")
source("covidmap/utils.r")
source("epimap/epimap.r")

rstan_options(auto_write = FALSE)

# Function defining the defualt arguments for the covidmap
# stage1 computations
covidmap_stage1_options = function(
  gp_time_scale        = 28.0, # units of 1 day
  gp_time_decay_scale  = .1,
  fixed_gp_time_length_scale = -1.0,

  first_day_modelled = NULL,
  last_day_modelled  = NULL,
  weeks_modelled     = NULL,
  days_ignored       = NULL,
  days_per_step      = 7,
  days_predicted     = 2,
  num_steps_forecasted = 3,

  num_samples        = 20,
  num_iterations     = 5000,
  num_chains         = 1,

  data_directory     = "data/",
  results_directory  = "fits/test",
  produce_plots      = FALSE,
  approximation      = "singlearea",
  area_index         = 0,

  limit_area         = NULL,
  limit_radius       = NULL,

  regions_as_areas_stage1   = FALSE,

  Adp                = 1.57,
  Aip                = 2.29,

  Bdp                = 0.65,
  Bip                = 0.36,

  num_bootstrap      = 10,

  bootstrap_merge      = FALSE

) {
  as.list(environment())
}
#' Produce plots to instect the results of the stage 1 model 
#' 
#' @param area index of the area to plot
#' @param Count The count data inference is done from
#' @param fit The stanfit from the stage 1 model
#' @param directory The directory to save the plots in 
#' @param Tcond The number of days conditioned on
#' @param Tstep The days per constant R step
#' @param Nstep The number of model steps
#' @param Nsample The number of samples take from the model
covidmap_stage1_plots = function(area,Count,fit,directory,Tcond,Tstep,Nstep,Nsample) {

  dir.create(directory, recursive=TRUE) 

  Tcur = Tcond + Tstep*Nstep
  ####################################################################
  # pairs plot
  pdf(paste(directory,"/pairs-", area, ".pdf", sep = ""), width = 9, height = 9)
  pairs(fit, pars = c("mu", "sigma", "alpha", "phi_latent", "phi_observed"))
  dev.off()

  Crecon_sample <- extract(fit, pars = "Crecon", permuted = FALSE)
  Clatent_sample <- extract(fit, pars = "Xt", permuted = FALSE)
  numsamples = dim(Crecon_sample)[1]
  skip <- numsamples / Nsample
  Clatent_sample <- Clatent_sample[seq(from=skip, by = skip, length.out = Nsample), , ]
  Crecon_sample <- Crecon_sample[seq(from=skip, by = skip, length.out = Nsample), , ]
  dim(Clatent_sample) <- c(Nsample, Tcur)
  Clatent_sample <- t(Clatent_sample)
  Crecon_sample <- t(Crecon_sample)
  dim(Clatent_sample) <- c(Tcur, Nsample)
  dim(Crecon_sample) <- c(Tcur, Nsample)
 
  pdf(paste(directory,"/Clatent-", area, ".pdf", sep = ""), width = 9, height = 9)
  par(mfrow = c(5, 2))
  par(oma = c(0, 0, 0, 0))
  par(mar = c(1, 1, 1, 1))
  ClatentCI <- summary(fit, pars = "Xt_proj", probs = c(0.025, 0.25, 0.5, 0.75, 0.975))$summary
  ind <- (Tcond + 1):Tcur
  ClatentCI <- ClatentCI[, c("2.5%", "50%", "97.5%")]


  for (i in 1:Nsample) {
    plot(t(Count), pch = 20, ylim = c(0, max(Count[, ind])))
    for (j in 1:3) {
      lines(ind, ClatentCI[ind, j])
    }
    points(ind, Clatent_sample[ind, i], col = "red", pch = 20)
  }
  dev.off()

  pdf(paste(directory,"/Crecon-", area, ".pdf", sep = ""), width = 9, height = 9)
  par(mfrow = c(5, 2))
  par(oma = c(0, 0, 0, 0))
  par(mar = c(1, 1, 1, 1))
  CreconCI <- summary(fit, pars = "Crecon", probs = c(0.025, 0.25, 0.5, 0.75, 0.975))$summary
  ind <- (Tcur - Nstep * Tstep + 1):Tcur
  CreconCI <- CreconCI[, c("2.5%", "50%", "97.5%")]
 
  for (i in 1:Nsample) {
    plot(t(Count), pch = 20, ylim = c(0, max(Count[ind])))
    for (j in 1:3) {
      lines(ind, CreconCI[ind, j])
    }
    points(ind, Crecon_sample[ind, i], col = "red", pch = "x")
  }
  dev.off()
}

#' Run a single area under the singlearea approximation
#' 
#' @param area_index The index of the area to run
#' @param opt A set options to run the model for, based on covidmap_stage1_options
covidmap_stage1_run = function(area_index = 0, opt = covidmap_stage1_options()) {
  if (area_index==0) {
    stop("Area index 0.")
  }
  env = covidmap_read_data(environment())

  numiters = opt$iterations


  # work out days to be modelled
  list[Nstep, Tstep, Tcond, Tlik, Tcur, Tignore] = process_dates_modelled(
    dates, 
    first_day_modelled = opt$first_day_modelled,
    last_day_modelled  = opt$last_day_modelled,
    days_ignored       = opt$days_ignored,
    weeks_modelled     = opt$weeks_modelled,
    days_per_step      = opt$days_per_step
  )

  area = areas[area_index]
  Count <- AllCount[area,]
  message("Area = ",area)

  Nsample <- opt$num_samples

  # Case only reported a few days after testing, 
  # no result delay truncation
  Trdp <- 1
  resultdelaydecay = .5
  resultdelaystrength = Trdp
  resultdelayalpha = resultdelaystrength * resultdelaydecay^(1:Trdp)

  # fit model using epimap
  start_time <- Sys.time()
  fit = epimap_singlearea(
    iter = opt$num_iterations,
    chains = opt$num_chains,
    cases = Count,
    Tcond = Tcond,
    Nstep = Nstep,
    Nproj = opt$num_steps_forecasted,
    Tstep = Tstep,
    Tpred = opt$days_predicted,
    infprofile = infprofile,
    testdelayprofile = testdelayprofile,
    resultdelayalpha = resultdelayalpha,
    gp_time_scale = opt$gp_time_scale,
    gp_time_decay_scale = opt$gp_time_decay_scale,
    fixed_gp_time_length_scale = opt$fixed_gp_time_length_scale
  )
  end_time <- Sys.time()
  
  # save results
  dir.create(paste(opt$results_directory,'/singlearea/stanfits',sep=''), recursive=TRUE) 
  saveRDS(fit, paste(opt$results_directory, '/singlearea/stanfits/',area,'.rds',sep=''))

  message(paste(
    "Time to run:",
    round(difftime(end_time, start_time, units="hours"),3),
    "hours")
  )

  # optionally prodice plots
  if (opt$produce_plots) {
    covidmap_stage1_plots(
      area,Count,fit,paste(opt$results_directory,"/singlearea/pdfs",sep=""),
      Tcond,Tstep,Nstep,Nsample
    )
  }
  return(invisible(fit))
}


covidmap_stage1_combine = function(opt = covidmap_stage1_options()) {
  covidmap_read_data(environment())

  # work out days to be modelled
  list[Nstep, Tstep, Tcond, Tlik, Tcur, Tignore] = process_dates_modelled(
    dates, 
    first_day_modelled = opt$first_day_modelled,
    last_day_modelled  = opt$last_day_modelled,
    days_ignored       = opt$days_ignored,
    weeks_modelled     = opt$weeks_modelled,
    days_per_step      = opt$days_per_step
  )

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

  Nproj = opt$num_steps_forecasted
  Tproj = Nproj * Tstep
  Tpred = opt$days_predicted
  provenance <- c(rep('inferred',Tlik),rep('projected',Tproj))
  days_likelihood = dates[(Tcond+1):Tcur]
  days_all <- c(days_likelihood,seq(days_likelihood[Tlik]+1,by=1,length.out=Tproj))
  message("Nstep = ",Nstep)
  message("Nproj = ",Nproj)

  Count <- AllCount[, 1:Tcur]

  numiters <- opt$num_iterations 
  Nsample <- opt$num_samples

  # Initialize arrays
  Clatent_sample <- array(0, c(N, Tcur+Tproj, Nsample))
  Clatent_mean <- array(0, c(N, Tcur+Tproj))
  Clatent_median <- array(0, c(N, Tcur+Tproj))
  Crecon_sample <- array(0, c(N, Tcur, Nsample))
  Crecon_median <- array(0, c(N, Tcur))
  Clatent_mean[, 1:Tcond] <- as.matrix(Count[, 1:Tcond])
  Clatent_median[, 1:Tcond] <- as.matrix(Count[, 1:Tcond])
  for (i in 1:Nsample) {
    Clatent_sample[, 1:Tcond, i] <- as.matrix(Count[, 1:Tcond])
    Crecon_sample[, 1:Tcond, i] <- as.matrix(Count[, 1:Tcond])
  }
  percentiles = c(.025,.1,.2,.25,.3,.4,.5,.6,.7,.75,.8,.9,.975)
  str_percentiles = c("2.5%","10%","20%","25%","30%","40%","50%","60%","70%","75%","80%","90%","97.5%")
  num_percentiles = length(percentiles)
  num_samples = opt$num_iterations
  Rt_percentiles = array(0.0, c(N, Nstep+Nproj, num_percentiles))
  Rt_samples = array(0.0, c(N, Nstep+Nproj, num_samples))
  logpred = array(0.0, c(N, Tpred))

  Bpred = array(0.0, c(N, Nstep*Tstep, num_percentiles))
  Bproj = array(0.0, c(N, Nproj*Tstep, num_percentiles))

  Cpred = array(0.0, c(N, Nstep*Tstep, num_percentiles))
  Cproj = array(0.0, c(N, Nproj*Tstep, num_percentiles))

  Xpred = array(0.0, c(N, Nstep*Tstep, num_percentiles))
  Xproj = array(0.0, c(N, Nproj*Tstep, num_percentiles))
  

  # Loop over areas, loading area RDS files and filling the arrays
  for (area_index in 1:N) {
    area <- areas[area_index]
    print(area)
  
    fit <- readRDS(paste(opt$results_directory, '/singlearea/stanfits/',area,'.rds',sep=''))
  
    skip <- numiters / 2 / Nsample
    ####################################################################
    Clatent_s <- extract(fit, pars = "Xt_proj", permuted = FALSE)
    Clatent_s <- Clatent_s[seq(from=skip, by = skip, length.out = Nsample), , ]
    dim(Clatent_s) <- c(Nsample, Tcur+Tproj)
    Clatent_s <- t(Clatent_s)
    dim(Clatent_s) <- c(1, Tcur+Tproj, Nsample)
    Clatent_sample[area_index, , ] <- Clatent_s
    Clatent_m <- summary(fit, pars = "Xt_proj", probs = c(0.5))$summary
    cmean <- t(as.matrix(Clatent_m[, "mean"]))
    Clatent_mean[area_index, ] <- cmean
    cmedian <- t(as.matrix(Clatent_m[, "50%"]))
    Clatent_median[area_index, ] <- cmedian

    ####################################################################
    Crecon_s <- extract(fit, pars = "Crecon", permuted = FALSE)
    Crecon_s <- Crecon_s[seq(from=skip, by = skip, length.out = Nsample), , ]
    Crecon_s <- t(Crecon_s)
    dim(Crecon_s) <- c(1, Tcur, Nsample)
    Crecon_sample[area_index, , ] <- Crecon_s
    Crecon_m <- summary(fit, pars = "Crecon", probs = c(0.5))$summary
    Crecon_m <- t(as.matrix(Crecon_m[, "50%"]))
    Crecon_median[area_index, ] <- round(Crecon_m)
  
    ####################################################################
    area_rt = summary(fit, pars = "Rt", probs=percentiles)$summary
    # message("area_rt dims = ",dim(area_rt))
    for (p in 1:num_percentiles) {
      Rt_percentiles[area_index,,p] = area_rt[1:(Nstep+Nproj),str_percentiles[p]]
    }
    samples = extract(fit, pars="Rt", permuted = FALSE)
    samples = samples[,1,]
    samples = aperm(samples, c(2,1))
    Rt_samples[area_index,,] = samples

    s <- summary(fit, pars="Ppred", probs=c(0.025, .5, .975))$summary
    Ppred <- s[,"mean"]
    logpred[area_index,] = log(Ppred)

    s <- summary(fit, pars="Bpred", probs=percentiles)$summary
    Bpred[area_index,,] = s[,str_percentiles]

    s <- summary(fit, pars="Bproj", probs=percentiles)$summary
    Bproj[area_index,,] = s[,str_percentiles]

    s <- summary(fit, pars="Cpred", probs=percentiles)$summary
    Cpred[area_index,,] = s[,str_percentiles]

    s <- summary(fit, pars="Cproj", probs=percentiles)$summary
    Cproj[area_index,,] = s[,str_percentiles]

    s <- summary(fit, pars="Xpred", probs=percentiles)$summary
    Xpred[area_index,,] = s[,str_percentiles]

    s <- summary(fit, pars="Xproj", probs=percentiles)$summary
    Xproj[area_index,,] = s[,str_percentiles]
 
  }
  
  days <- colnames(Count)
  days_proj <- c(days,as.character(seq(days_likelihood[Tlik]+1,by=1,length.out=Tproj),format='%Y-%m-%d'))

  rownames(Clatent_mean) <- quoted_areas
  colnames(Clatent_mean) <- days_proj
  write.csv(Clatent_mean, paste(opt$results_directory, "/singlearea/Clatent_mean.csv", sep=""), quote = FALSE)
  rownames(Clatent_median) <- quoted_areas
  colnames(Clatent_median) <- days_proj
  write.csv(Clatent_median, paste(opt$results_directory, "/singlearea/Clatent_median.csv", sep=""), quote = FALSE)
  rownames(Crecon_median) <- quoted_areas
  colnames(Crecon_median) <- days
  write.csv(Crecon_median, paste(opt$results_directory, "/singlearea/Crecon_median.csv", sep=""), quote = FALSE)
  for (i in 1:Nsample) {
    cc <- Clatent_sample[, , i]
    rownames(cc) <- quoted_areas
    colnames(cc) <- days_proj
    write.csv(cc, paste(opt$results_directory, "/singlearea/Clatent_sample", i, ".csv", sep = ""), quote = FALSE)
    cc <- Crecon_sample[, , i]
    rownames(cc) <- quoted_areas
    colnames(cc) <- days
    write.csv(cc, paste(opt$results_directory, "/singlearea/Crecon_sample", i, ".csv", sep = ""), quote = FALSE)
  }

  # days_modelled = days[seq(from=Tcond+1, by=Tstep, to=Tcond+Tlik)]
  # Rt_percentiles = format(round(Rt_percentiles,2),nsmall=2)
  # dimnames(Rt_percentiles) = list(area=quoted_areas,date=days_modelled,str_percentiles)
  # df = adply(Rt_percentiles,c(2,1))
  # df = df[,c(2,1,3:(num_percentiles+2))]

  Rt = Rt_percentiles[,rep(c(1:(Nstep+Nproj)),each=Tstep),]
  Rt = aperm(Rt, c(2,1,3))
  dim(Rt) <- c(N*(Nstep+Nproj)*Tstep, num_percentiles)
  Rt <- area_date_dataframe(
    quoted_areas,
    days_all,
    provenance,
    format(round(Rt,2),nsmall=2),
    c("Rt_2_5","Rt_10","Rt_20","Rt_25","Rt_30","Rt_40","Rt_50",
      "Rt_60","Rt_70","Rt_75","Rt_80","Rt_90","Rt_97_5")
  )
  write.csv(Rt, paste(opt$results_directory, "/singlearea/Rt.csv", sep=""), quote=FALSE, row.names=FALSE)

  thresholds = c(.8, .9, 1.0, 1.1, 1.2, 1.5, 2.0)
  numthresholds = length(thresholds)
  Pexceedance = array(0.0, dim=c(Nstep+Nproj,N,numthresholds))
  for (k in 1:(Nstep+Nproj)) {
    for (i in 1:N) {
      for (x in 1:numthresholds) {
        Pexceedance[k,i,x] = mean(Rt_samples[i,k,]>thresholds[x])
      }
    }
  }
  Pexceedance = Pexceedance[c(1:(Nstep+Nproj)),,]
  Pexceedance <- Pexceedance[sapply(1:(Nstep+Nproj),function(k)rep(k,Tstep)),,]
  dim(Pexceedance) <- c(Tstep*(Nstep+Nproj)*N,numthresholds)
  Pexceedance <- area_date_dataframe(
      quoted_areas, 
      days_all,
      provenance,
      format(round(Pexceedance,2),nsmall=2),
      c("P_08","P_09","P_10","P_11","P_12","P_15","P_20")
  )
  write.csv(Pexceedance, paste(opt$results_directory, "/singlearea/Pexceed.csv", sep=""), quote=FALSE, row.names=FALSE)

  logpred = data.frame(area = areas, logpred = logpred, provenance=rep('inferred', length(areas)))
  for (i in 1:Tpred)
      colnames(logpred)[i+1] <- sprintf('logpred_day%d',i)
  write.csv(logpred, paste(opt$results_directory, "/singlearea/logpred.csv", sep=""), quote=FALSE, row.names=FALSE)


  Bweekly = array(0.0, c(N, (Nstep + Nproj), Tstep))

  actuals <- as.matrix(AllCount[,(Tcond+1):(Tcond+Tlik)])
  dim(actuals) <- c(N,Tstep,Nstep)
  actuals <- aperm(actuals, c(1,3,2))
  # preds = Cpred[,,3]
  # dim(preds) <- c(N, Nstep, Tstep)
  Bweekly[,1:Nstep,] = actuals

  projs = Bproj[,,3]
  dim(projs) <- c(N, Tstep, Nproj)
  projs = aperm(projs,c(1,3,2))
  Bweekly[,(Nstep+1):(Nstep+Nproj),] = projs

  Bweekly = apply(Bweekly, c(1,2), sum) 
  Bweekly = Bweekly[,sapply(1:(Nstep+Nproj),function(k)rep(k,Tstep))]
  Bweekly = aperm(Bweekly, c(2,1))
  dim(Bweekly) <- c(N*(Nstep+Nproj)*Tstep)

  Bweekly <- area_date_dataframe(
    quoted_areas,
    days_all,
    provenance,
    Bweekly,
    c("B_weekly")
  )
  write.csv(Bweekly, paste(opt$results_directory, "/singlearea/Bweekly.csv", sep=""), quote=FALSE, row.names=FALSE)

  
  Cweekly = array(0.0, c(N, (Nstep + Nproj), Tstep))

  actuals <- as.matrix(AllCount[,(Tcond+1):(Tcond+Tlik)])
  dim(actuals) <- c(N,Tstep,Nstep)
  actuals <- aperm(actuals, c(1,3,2))
  # preds = Cpred[,,3]
  # dim(preds) <- c(N, Nstep, Tstep)
  Cweekly[,1:Nstep,] = actuals

  projs = Cproj[,,3]
  dim(projs) <- c(N, Tstep, Nproj)
  projs = aperm(projs,c(1,3,2))
  Cweekly[,(Nstep+1):(Nstep+Nproj),] = projs

  Cweekly = apply(Cweekly, c(1,2), sum) 
  Cweekly = Cweekly[,sapply(1:(Nstep+Nproj),function(k)rep(k,Tstep))]
  Cweekly = aperm(Cweekly, c(2,1))
  dim(Cweekly) <- c(N*(Nstep+Nproj)*Tstep)

  Cweekly <- area_date_dataframe(
    quoted_areas,
    days_all,
    provenance,
    Cweekly,
    c("C_weekly")
  )
  write.csv(Cweekly, paste(opt$results_directory, "/singlearea/Cweekly.csv", sep=""), quote=FALSE, row.names=FALSE)


  Bpred = aperm(Bpred, c(2,1,3))
  dim(Bpred) <- c(N*Nstep*Tstep, num_percentiles)
  Bpred <- area_date_dataframe(
    quoted_areas,
    days_likelihood,
    rep('inferred',Nstep*Tstep),
    Bpred,
    c("C_025","C_10","C_20","C_25","C_30","C_40","C_50",
        "C_60", "C_70","C_75","C_80","C_90","C_975")
  )
  write.csv(Bpred, paste(opt$results_directory, "/singlearea/Bpred.csv", sep=""), quote=FALSE, row.names=FALSE)


  Bproj = aperm(Bproj, c(2,1,3))
  dim(Bproj) <- c(N*Nproj*Tstep, num_percentiles)
  Bproj <- area_date_dataframe(
    quoted_areas,
    seq(days_likelihood[Tlik]+1,by=1,length.out=Tproj),
    rep('projected',Tproj),
    Bproj,
    c("C_025","C_10","C_20","C_25","C_30","C_40","C_50",
        "C_60", "C_70","C_75","C_80","C_90","C_975")
  )
  write.csv(Bproj, paste(opt$results_directory, "/singlearea/Bproj.csv", sep=""), quote=FALSE, row.names=FALSE)


  Cpred = aperm(Cpred, c(2,1,3))
  dim(Cpred) <- c(N*Nstep*Tstep, num_percentiles)
  Cpred <- area_date_dataframe(
    quoted_areas,
    days_likelihood,
    rep('inferred',Nstep*Tstep),
    Cpred,
    c("C_025","C_10","C_20","C_25","C_30","C_40","C_50",
        "C_60", "C_70","C_75","C_80","C_90","C_975")
  )
  write.csv(Cpred, paste(opt$results_directory, "/singlearea/Cpred.csv", sep=""), quote=FALSE, row.names=FALSE)


  Cproj = aperm(Cproj, c(2,1,3))
  dim(Cproj) <- c(N*Nproj*Tstep, num_percentiles)
  Cproj <- area_date_dataframe(
    quoted_areas,
    seq(days_likelihood[Tlik]+1,by=1,length.out=Tproj),
    rep('projected',Tproj),
    Cproj,
    c("C_025","C_10","C_20","C_25","C_30","C_40","C_50",
        "C_60", "C_70","C_75","C_80","C_90","C_975")
  )
  write.csv(Cproj, paste(opt$results_directory, "/singlearea/Cproj.csv", sep=""), quote=FALSE, row.names=FALSE)

  Xpred = aperm(Xpred, c(2,1,3))
  dim(Xpred) <- c(N*Nstep*Tstep, num_percentiles)
  Xpred <- area_date_dataframe(
    quoted_areas,
    days_likelihood,
    rep('inferred',Nstep*Tstep),
    Xpred,
    c("X_025","X_10","X_20","X_25","X_30","X_40","X_50",
        "X_60", "X_70","X_75","X_80","X_90","X_975")
  )
  write.csv(Xpred, paste(opt$results_directory, "/singlearea/Xpred.csv", sep=""), quote=FALSE, row.names=FALSE)


  Xproj = aperm(Xproj, c(2,1,3))
  dim(Xproj) <- c(N*Nproj*Tstep, num_percentiles)
  Xproj <- area_date_dataframe(
    quoted_areas,
    seq(days_likelihood[Tlik]+1,by=1,length.out=Tproj),
    rep('projected',Tproj),
    Xproj,
    c("X_025","X_10","X_20","X_25","X_30","X_40","X_50",
        "X_60", "X_70","X_75","X_80","X_90","X_975")
  )
  write.csv(Xproj, paste(opt$results_directory, "/singlearea/Xproj.csv", sep=""), quote=FALSE, row.names=FALSE)


  list(
    Rt = Rt_percentiles,
    Clatent_sample = Clatent_sample,
    Clatent_mean = Clatent_mean,
    Crecon_sample = Crecon_sample,
    Crecon_median = Crecon_median
  )
}


covidmap_stage1_bootstrap_combine = function(opt = covidmap_stage1_options()) {
  covidmap_read_data(environment())

  # work out days to be modelled
  list[Nstep, Tstep, Tcond, Tlik, Tcur, Tignore] = process_dates_modelled(
    dates, 
    first_day_modelled = opt$first_day_modelled,
    last_day_modelled  = opt$last_day_modelled,
    days_ignored       = opt$days_ignored,
    weeks_modelled     = opt$weeks_modelled,
    days_per_step      = opt$days_per_step
  )

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

  Nproj = opt$num_steps_forecasted
  Tproj = Nproj * Tstep
  Tpred = opt$days_predicted
  provenance <- c(rep('inferred',Tlik),rep('projected',Tproj))
  days_likelihood = dates[(Tcond+1):Tcur]
  days_all <- c(days_likelihood,seq(days_likelihood[Tlik]+1,by=1,length.out=Tproj))
  message("Nstep = ",Nstep)
  message("Nproj = ",Nproj)

  Count <- AllCount[, 1:Tcur]

  numiters <- opt$num_iterations 
  Nsample <- opt$num_samples
  Nbootstrap <- opt$num_bootstrap

  percentiles = c(.025,.1,.2,.25,.3,.4,.5,.6,.7,.75,.8,.9,.975)
  str_percentiles = c("2.5%","10%","20%","25%","30%","40%","50%","60%","70%","75%","80%","90%","97.5%")
  num_percentiles = length(percentiles)
  num_samples = opt$num_iterations

  Rt_samples = array(0.0, c(Nbootstrap, Nstep+Nproj, num_samples))
  Cpred_samples = array(0.0, c(Nbootstrap, Nstep*Tstep, num_samples))
  Cproj_samples = array(0.0, c(Nbootstrap, Nproj*Tstep, num_samples))
  Rt_percentiles = array(0.0, c(N, Nstep+Nproj, num_percentiles))
  Cproj = array(0.0, c(N, Nproj*Tstep, num_percentiles))
  Cpred = array(0.0, c(N, Nstep*Tstep, num_percentiles))

  # Loop over areas, loading area RDS files and filling the arrays
  for (area_index in 1:2) {
    area <- areas[area_index]
    print(area)
    for (bootstrap_id in 1:Nbootstrap) {
      fit <- readRDS(paste(opt$results_directory, '/bootstrap_', bootstrap_id, '/singlearea/stanfits/',area,'.rds',sep=''))

      samples = extract(fit, pars="Rt", permuted = FALSE)
      samples = samples[,1,]
      samples = aperm(samples, c(2,1))
      Rt_samples[bootstrap_id,,] = samples

      samples = extract(fit, pars="Cpred", permuted = FALSE)
      samples = samples[,1,]
      samples = aperm(samples, c(2,1))
      Cpred_samples[bootstrap_id,,] = samples

      samples = extract(fit, pars="Cproj", permuted = FALSE)
      samples = samples[,1,]
      samples = aperm(samples, c(2,1))
      Cproj_samples[bootstrap_id,,] = samples
    }
    Rt_percentiles[area_index,,] = aperm(apply(Rt_samples, 2, quantile, probs=percentiles), c(2,1))
    Cpred[area_index,,] = aperm(apply(Cpred_samples, 2, quantile, probs=percentiles), c(2,1))
    Cproj[area_index,,] = aperm(apply(Cproj_samples, 2, quantile, probs=percentiles), c(2,1))
  }
  
  days <- colnames(Count)
  days_proj <- c(days,as.character(seq(days_likelihood[Tlik]+1,by=1,length.out=Tproj),format='%Y-%m-%d'))

  Rt = Rt_percentiles[,rep(c(1:(Nstep+Nproj)),each=Tstep),]
  Rt = aperm(Rt, c(2,1,3))
  dim(Rt) <- c(N*(Nstep+Nproj)*Tstep, num_percentiles)
  Rt <- area_date_dataframe(
    quoted_areas,
    days_all,
    provenance,
    format(round(Rt,2),nsmall=2),
    c("Rt_2_5","Rt_10","Rt_20","Rt_25","Rt_30","Rt_40","Rt_50",
      "Rt_60","Rt_70","Rt_75","Rt_80","Rt_90","Rt_97_5")
  )
  write.csv(Rt, paste(opt$results_directory, "/singlearea/Rt.csv", sep=""), quote=FALSE, row.names=FALSE)
  message("Saved Rt")

  Cpred = aperm(Cpred, c(2, 1, 3))
  dim(Cpred) <- c(N*Nstep*Tstep, num_percentiles)
  Cpred <- area_date_dataframe(
    quoted_areas,
    days_likelihood,
    rep('inferred',Nstep*Tstep),
    Cpred,
    c("C_025","C_10","C_20","C_25","C_30","C_40","C_50",
        "C_60", "C_70","C_75","C_80","C_90","C_975")
  )
  write.csv(Cpred, paste(opt$results_directory, "/singlearea/Cpred.csv", sep=""), quote=FALSE, row.names=FALSE)
  message("Saved Cpred")

  Cproj = aperm(Cproj, c(2, 1, 3))
  dim(Cproj) <- c(N*Nproj*Tstep, num_percentiles)
  Cproj <- area_date_dataframe(
    quoted_areas,
    seq(days_likelihood[Tlik]+1,by=1,length.out=Tproj),
    rep('projected',Tproj),
    Cproj,
    c("C_025","C_10","C_20","C_25","C_30","C_40","C_50",
        "C_60", "C_70","C_75","C_80","C_90","C_975")
  )
  write.csv(Cproj, paste(opt$results_directory, "/singlearea/Cproj.csv", sep=""), quote=FALSE, row.names=FALSE)
  message("Saved Cproj")
}

covidmap_stage1_cmdline_options = function(opt = covidmap_stage1_options()) {
  list(
    make_option(
      c("--area_index"),
      type="integer",
      default=0,
      help="Area index (required argument)."
    ),
    make_option(
      c("--num_samples"),
      type="integer",
      default=opt$num_samples,
      help="Number of samples to output."
    ),
    make_option(
      c("--num_iterations"),
      type="integer",
      default=opt$num_iterations,
      help=paste("Number of MCMC iterations, default",opt$num_iterations)
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
      help=paste("Number of recent days ignored, default",opt$days_ignored)
    ),
    make_option(
      c("--days_per_step"),
      type="integer",
      default=opt$days_per_step,
      help=paste("Days per modelling step, default",opt$days_per_step)
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
      help=paste("Number of steps to forecast, default",opt$num_steps_forecasted)
    ),
    make_option(
      c("--gp_time_scale"),
      type="double",
      default=opt$gp_time_scale,
      help=paste("GP time scale, default",opt$gp_time_scale)
    ),
    make_option(
      c("--gp_time_decay_scale"),
      type="double",
      default=opt$gp_time_decay_scale,
      help=paste("GP time decay scale, default",opt$gp_time_decay_scale)
    ),
    make_option(
      c("--fixed_gp_time_length_scale"),
      type="double",
      default=opt$fixed_gp_time_length_scale,
      help=paste("Fixed GP time length scale, default",
                 opt$fixed_gp_time_length_scale)
    ),

    make_option(
      c("--results_directory"), 
      type="character", 
      default=opt$results_directory, 
      help=paste("Directory to put cleaned results in, default ",opt$results_directory)
    ),

    make_option(
      c("--data_directory"), 
      type="character", 
      default=opt$data_directory, 
      help=paste("Directory to get data from, default ",opt$data_directory)
    ),

    make_option(
      c("--produce_plots"), 
      type="logical", 
      default=opt$produce_plots, 
      help=paste("Whether to produce plots; default",opt$produce_plots)
    ),

    make_option(
      c("--regions_as_areas_stage1"), 
      type="logical", 
      default=opt$regions_as_areas_stage1, 
      help=paste("Whether to run stage 1 on regions instead of areas; default", opt$regions_as_areas_stage1)
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
    ),

    make_option(
      c("--Adp"), 
      type="double", 
      default=opt$Adp, 
      help=paste("Shape parameter of test delay profile gamma dist; default", opt$Adp)
    ),
    make_option(
      c("--Aip"), 
      type="double", 
      default=opt$Aip, 
      help=paste("Shape parameter of generation interval gamma dist; default", opt$Aip)
    ),

    make_option(
      c("--Bdp"), 
      type="double", 
      default=opt$Bdp, 
      help=paste("Scale parameter of test delay profile gamma dist; default", opt$Bdp)
    ),
    make_option(
      c("--Bip"), 
      type="double", 
      default=opt$Bip, 
      help=paste("Scale parameter of generation interval gamma dist; default", opt$Bip)
    ),
    make_option(
      c("--num_bootstrap"),
      type="integer",
      default=1,
      help="Number of bootstrap samples."
    )
  )
}

covidmap_stage1_get_cmdline_options = function(opt=covidmap_stage1_options()) {
  cmdline_opt = covidmap_stage1_cmdline_options(opt)
  opt_parser = OptionParser(option_list=cmdline_opt)
  parsed_opt = parse_args(opt_parser)
  for (o in names(parsed_opt)) {
    opt[o] = parsed_opt[o]
  }
  opt
}


