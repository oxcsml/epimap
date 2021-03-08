library(EpiEstim)
library(optparse)
library(gsubfn)
library(plyr)
library(data.table) 
source("covidmap/read_data.r")
source("covidmap/utils.r")


epiestim_options = function(
  first_day_modelled   = NULL,
  last_day_modelled    = NULL,
  weeks_modelled       = 15,
  days_ignored         = 7,
  days_per_step        = 7,
  days_predicted       = 0,
  num_steps_forecasted = 0,

  num_samples        = 20,
  num_iterations     = 3000,
  num_chains         = 1,

  data_directory     = "data/",
  results_directory  = "fits/test",
  produce_plots      = FALSE,
  approximation      = "epiestim",
  area_index         = 0
) {
  as.list(environment())
}


epiestim_run = function(area_index = 0, opt = epiestim_options()) {
  if (area_index==0) {
    stop("Area index 0.")
  }
  env = covidmap_read_data(environment())

  numiters = opt$iterations

  list[Nstep, Tstep, Tcond, Tlik, Tcur, Tignore] = process_dates_modelled(
    dates, 
    first_day_modelled = opt$first_day_modelled,
    last_day_modelled  = opt$last_day_modelled,
    days_ignored       = opt$days_ignored,
    weeks_modelled     = opt$weeks_modelled,
    days_per_step      = opt$days_per_step
  )
  message("Nproj = ", opt$num_steps_forecasted)

  area = areas[area_index]
  message("Area = ",area)

  Nsample <- opt$num_samples

  # Count = unlist(AllCount[area,(Tcond):(Tcond+(Tstep*Nstep))], use.names = FALSE)
  Count = AllCount[area,(Tcond):(Tcond+(Tstep*Nstep))]
  dates = colnames(Count)
  Count = transpose(Count)
  starts = seq(from=2,by=Tstep,length.out=Nstep) 
  ends = seq(from=Tstep+1,by=Tstep,length.out=Nstep)

#   starts = seq(from=2,by=1,length.out=(Nstep*Tstep)-Tstep-1) 
#   ends = seq(from=Tstep+1,by=1,length.out=(Nstep*Tstep)-Tstep-1)

  modified_infprofile = replicate(length(infprofile)+1, 0.0)
  modified_infprofile[2:(length(infprofile)+1)] = infprofile
  modified_infprofile = unlist(modified_infprofile)

  config = make_config(
    t_start = starts, 
    t_end = ends, 
    si_distr=modified_infprofile,
    mean_prior=1.0
  )

  fit = estimate_R(
    Count, 
    method="non_parametric_si",
    config = config
  )
  wt = wallinga_teunis(
    Count, 
    method="non_parametric_si",
    config = config
  )

  fit$dates = dates
  dir.create(paste(opt$results_directory,'/epiestim/fits',sep=''), recursive=TRUE) 
  saveRDS(fit, paste(opt$results_directory, '/epiestim/fits/',area,'.rds',sep=''))

}


epiestim_combine = function(opt = epiestim_options()) {
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

  # Initialize dummy arrays
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
  logpred = array(0.0, c(N, Tpred))

  C_percentiles = c(.025,.25,.5,.75,.975)
  C_str_percentiles = c("2.5%","25%","50%","75%","97.5%")
  num_C_percentiles = length(C_percentiles)
  Cpred = array(0.0, c(N, Nstep*Tstep, num_C_percentiles))
  Cproj = array(0.0, c(N, Nproj*Tstep, num_C_percentiles))
  

  # Initialise actual arrays
  # percentiles = c(.025,.1,.2,.25,.3,.4,.5,.6,.7,.75,.8,.9,.975)
  # str_percentiles = c("2.5%","10%","20%","25%","30%","40%","50%","60%","70%","75%","80%","90%","97.5%")
  percentiles = c(.025,.05,.25,.5,.75,.095,.975)
  str_percentiles = c("2.5%","5%","25%","50%","75%","95%","97.5%")
  epiestim_percentiles = c("Quantile.0.025(R)","Quantile.0.05(R)","Quantile.0.25(R)","Median(R)","Quantile.0.75(R)","Quantile.0.95(R)","Quantile.0.975(R)")
  # pad_percentiles = c("10%", "20%", "30%", "40%", "60%", "70%", "80%", "90%")
  pad_percentiles = c()
  num_percentiles = length(percentiles)
  num_samples = opt$num_iterations
  Rt_percentiles = array(0.0, c(N, Nstep+Nproj, num_percentiles + length(pad_percentiles)))
  Rt_samples = array(0.0, c(N, Nstep+Nproj, num_samples))

  # Loop over areas, loading area RDS files and filling the arrays
  for (area_index in 1:N) {
    area <- areas[area_index]
    print(area)
  
    fit <- readRDS(paste(opt$results_directory, '/epiestim/fits/',area,'.rds',sep=''))
    for (p in 1:num_percentiles) {
      Rt_percentiles[area_index,1:Nstep, p] = fit$R[,epiestim_percentiles[p]]
      for (i in 1:Nproj) {
        Rt_percentiles[area_index,Nstep + i, p] = fit$R[Nstep,epiestim_percentiles[p]]
      }
    }
    
  }

  days <- colnames(Count)
  days_proj <- c(days,as.character(seq(days_likelihood[Tlik]+1,by=1,length.out=Tproj),format='%Y-%m-%d'))

  # Save real data
  Rt = Rt_percentiles[,rep(c(1:(Nstep+Nproj)),each=Tstep),]
  Rt = aperm(Rt, c(2,1,3))
  dim(Rt) <- c(N*(Nstep+Nproj)*Tstep, num_percentiles + length(pad_percentiles))
  Rt <- area_date_dataframe(
    quoted_areas,
    days_all,
    provenance,
    format(round(Rt,2),nsmall=2),
    c("Rt_2_5","Rt_5","Rt_25","Rt_50","Rt_75","Rt_95","Rt_97_5") #, "Rt_10", "Rt_20", "Rt_30", "Rt_40", "Rt_60", "Rt_70", "Rt_80", "Rt_90")
  )
  write.csv(Rt, paste(opt$results_directory, "/epiestim/Rt.csv", sep=""), quote=FALSE, row.names=FALSE)


  # Save dummy data
  thresholds = c(.8, .9, 1.0, 1.1, 1.2, 1.5, 2.0)
  numthresholds = length(thresholds)
  Pexceedance = array(0.0, dim=c(Nstep+Nproj,N,numthresholds))
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
  write.csv(Pexceedance, paste(opt$results_directory, "/epiestim/Pexceed.csv", sep=""), quote=FALSE, row.names=FALSE)


  logpred = data.frame(area = areas, logpred = logpred, provenance=rep('inferred', length(areas)))
  for (i in 1:Tpred)
      colnames(logpred)[i+1] <- sprintf('logpred_day%d',i)
  write.csv(logpred, paste(opt$results_directory, "/epiestim/logpred.csv", sep=""), quote=FALSE, row.names=FALSE)

  Cweekly = array(0.0, c(N, (Nstep + Nproj), Tstep))
  actuals <- as.matrix(AllCount[,(Tcond+1):(Tcond+Tlik)])
  dim(actuals) <- c(N,Tstep,Nstep)
  actuals <- aperm(actuals, c(1,3,2))
  # preds = Cpred[,,3]
  # dim(preds) <- c(N, Nstep, Tstep)
  Cweekly[,1:Nstep,] = actuals
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
  write.csv(Cweekly, paste(opt$results_directory, "/epiestim/Cweekly.csv", sep=""), quote=FALSE, row.names=FALSE)

  Cpred = aperm(Cpred, c(2,1,3))
  dim(Cpred) <- c(N*Nstep*Tstep, num_C_percentiles)
  Cpred <- area_date_dataframe(
    quoted_areas,
    days_likelihood,
    rep('inferred',Nstep*Tstep),
    Cpred,
    c("C_025","C_25","C_50","C_75","C_975")
  )
  write.csv(Cpred, paste(opt$results_directory, "/epiestim/Cpred.csv", sep=""), quote=FALSE, row.names=FALSE)

  Cproj = aperm(Cproj, c(2,1,3))
  dim(Cproj) <- c(N*Nproj*Tstep, num_C_percentiles)
  Cproj <- area_date_dataframe(
    quoted_areas,
    seq(days_likelihood[Tlik]+1,by=1,length.out=Tproj),
    rep('projected',Tproj),
    Cproj,
    c("C_025","C_25","C_50","C_75","C_975")
  )
  write.csv(Cproj, paste(opt$results_directory, "/epiestim/Cproj.csv", sep=""), quote=FALSE, row.names=FALSE)


}

epiestim_cmdline_options = function(opt = epiestim_options()) {
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

epiestim_get_cmdline_options = function(opt=epiestim_options()) {
  cmdline_opt = epiestim_cmdline_options(opt)
  opt_parser = OptionParser(option_list=cmdline_opt)
  parsed_opt = parse_args(opt_parser)
  for (o in names(parsed_opt)) {
    opt[o] = parsed_opt[o]
  }
  opt
}

