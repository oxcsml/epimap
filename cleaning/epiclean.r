library(rstan)
library(optparse)
library(gsubfn)
library(plyr)
source("dataprocessing/read_data.r")
source("mapping/utils.r")
rstan_options(auto_write = FALSE)

epiclean_options = function(
  gp_time_scale        = 14.0, # units of 1 day
  gp_time_decay_scale  = .1,
  fixed_gp_time_length_scale = -1.0,

  first_day_modelled = "2020-06-01",
  last_day_modelled  = NULL,
  weeks_modelled     = NULL,
  days_ignored         = 7,
  days_per_step      = 1,
  num_steps_forecasted = 3*7,

  num_samples        = 20,
  iterations         = 3000,

  data_directory     = "data/",
  clean_directory    = "fits/clean",
  produce_plots      = FALSE
) {
  as.list(environment())
}

epiclean = function(area_index = 0, opt = epiclean_options()) {
  env = new.env(parent=globalenv())
  env$area_index = area_index
  if (area_index==0) {
    stop("Area index 0.")
  }
  env$opt = opt
  Rmap_read_data(env)
  with(env,{

    numchains = 1
    numiters = opt$iterations

    options(mc.cores = min(numchains,parallel::detectCores()))

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

    # ------------------------------------------------------------------------ #
    #                             Main computation                             #
    # ------------------------------------------------------------------------ #

    Rmap_clean_data <- list(
      Tall = Tall,
      Tcond = Tcond,
      Tstep = Tstep, 
      Nstep = Nstep,
      Nproj = opt$num_steps_forecasted,
      Count = Count[area,],
      Tip = Tip,
      infprofile = infprofile,
      Ttdp = length(testdelayprofile),
      testdelayprofile = testdelayprofile,
      Trdp = Trdp,
      resultdelaydecay = resultdelaydecay,
      resultdelaystrength = resultdelaystrength,
      gp_time_scale        = opt$gp_time_scale,
      gp_time_decay_scale  = opt$gp_time_decay_scale,
      fixed_gp_time_length_scale = opt$fixed_gp_time_length_scale,
      mu_scale = 0.5,
      sigma_scale = 0.5,
      # alpha_scale = 1-exp(-Tstep/14),
      phi_latent_scale = 1e-6,
      phi_observed_scale = 5.0,
      outlier_prob_threshold = .95,
      outlier_count_threshold = 2,
      xi_scale = 0.01,
      reconstruct_infections = TRUE
    )
    
    init = list()
    init[[1]] = list(
      mu = 0.0,
      sigma = 0.01,
      alpha1 = .95,
      phi_latent = 1.0,
      phi_observed = 5.0,
      xi = .01
    )
    for (i in 1:Nstep) {
      init[[1]][[paste('Reta[',i,']',sep='')]] = 0.0
    }
    
    start_time <- Sys.time()
    
    fit <- stan(file = 'cleaning/stan_files/Rmap-clean.stan',
                data = Rmap_clean_data, 
                init = init,
                iter = numiters, 
                chains = numchains,
                control = list(adapt_delta = .9))
    
    end_time <- Sys.time()
    
    print("Time to run")
    print(end_time - start_time)
   
    dir.create(paste(opt$clean_directory,'/stanfits',sep=''), showWarnings = FALSE) 
    saveRDS(fit, paste(opt$clean_directory, '/stanfits/',area,'.rds',sep=''))
    
    #################################################################
    # Summary of fit
    print(summary(fit, 
      pars=c("mu","sigma","alpha","gp_time_length_scale","phi_latent","phi_observed","xi","Noutliers","meandelay","resultdelayprofile","Rt"),
      probs=c(0.5)
    )$summary)
  })
  env
}

epiclean_combine = function(opt = epiclean_options()) {
  env = new.env(parent=globalenv())
  env$opt = opt
  Rmap_read_data(env)
  with(env,{

    # work out days to be modelled
    list[Nstep, Tstep, Tcond, Tlik, Tcur, Tignore] = process_dates_modelled(
      dates, 
      first_day_modelled = opt$first_day_modelled,
      last_day_modelled  = opt$last_day_modelled,
      days_ignored       = opt$days_ignored,
      weeks_modelled     = opt$weeks_modelled,
      days_per_step      = opt$days_per_step
    )

    Count <- AllCount[, 1:Tcur]

    numiters <- opt$iterations 
    Nsample <- opt$num_samples

    # Initialize arrays
    Clatent_sample <- array(0, c(N, Tcur, Nsample))
    Clatent_mean <- array(0, c(N, Tcur))
    Crecon_sample <- array(0, c(N, Tcur, Nsample))
    Crecon_median <- array(0, c(N, Tcur))
    Clatent_mean[, 1:Tcond] <- as.matrix(Count[, 1:Tcond])
    for (i in 1:Nsample) {
      Clatent_sample[, 1:Tcond, i] <- as.matrix(Count[, 1:Tcond])
      Crecon_sample[, 1:Tcond, i] <- as.matrix(Count[, 1:Tcond])
    }
    percentiles = c(.025,.25,.5,.75,.975)
    str_percentiles = c("2.5%","25%","50%","75%","97.5%")
    num_percentiles = length(percentiles)
    Rt_percentiles = array(0, c(N, Tlik, num_percentiles))
    
    # Loop over areas, loading area RDS files and filling the arrays
    for (area_index in 1:N) {
      area <- areas[area_index]
      print(area)
    
      fit <- readRDS(paste(opt$clean_directory, '/stanfits/',area,'.rds',sep=''))
    
      skip <- numiters / 2 / Nsample
      ####################################################################
      Clatent_s <- extract(fit, pars = "Clatent", permuted = FALSE)
      Clatent_s <- Clatent_s[seq(skip, by = skip, length.out = Nsample), , ]
      dim(Clatent_s) <- c(Nsample, Tcur)
      Clatent_s <- t(Clatent_s)
      dim(Clatent_s) <- c(1, Tcur, Nsample)
      Clatent_sample[area_index, , ] <- Clatent_s
      Clatent_m <- summary(fit, pars = "Clatent", probs = c(0.5))$summary
      Clatent_m <- t(as.matrix(Clatent_m[, "mean"]))
      Clatent_mean[area_index, ] <- Clatent_m

      ####################################################################
      Crecon_s <- extract(fit, pars = "Crecon", permuted = FALSE)
      Crecon_s <- Crecon_s[seq(skip, by = skip, length.out = Nsample), , ]
      Crecon_s <- t(Crecon_s)
      dim(Crecon_s) <- c(1, Tcur, Nsample)
      Crecon_sample[area_index, , ] <- Crecon_s
      Crecon_m <- summary(fit, pars = "Crecon", probs = c(0.5))$summary
      Crecon_m <- t(as.matrix(Crecon_m[, "50%"]))
      Crecon_median[area_index, ] <- round(Crecon_m)

    
      ####################################################################
      area_rt = summary(fit, pars = "Rt", probs=percentiles)$summary
      for (p in 1:num_percentiles) {
        Rt_percentiles[area_index,,p] = area_rt[1:Tlik,str_percentiles[p]]
      }

      if (opt$produce_plots) {
        ####################################################################
        # pairs plot
        pdf(paste(opt$clean_directory,"/pdfs/pairs-", area, ".pdf", sep = ""), width = 9, height = 9)
        pairs(fit, pars = c("mu", "sigma", "alpha", "phi_latent", "phi_observed"))
        dev.off()
    
        pdf(paste(opt$clean_directory,"/pdfs/Clatent-", area, ".pdf", sep = ""), width = 9, height = 9)
        par(mfrow = c(5, 2))
        par(oma = c(0, 0, 0, 0))
        par(mar = c(1, 1, 1, 1))
        ClatentCI <- summary(fit, pars = "Clatent", probs = c(0.025, 0.25, 0.5, 0.75, 0.975))$summary
        ind <- (Tcond + 1):Tcur
        ClatentCI <- ClatentCI[, c("2.5%", "50%", "97.5%")]
        for (i in 1:Nsample) {
          plot(t(Count[area, ]), pch = 20, ylim = c(0, max(Count[area, ind])))
          for (j in 1:3) {
            lines(ind, ClatentCI[ind, j])
          }
          points(ind, Clatent_sample[area_index, ind, i], col = "red", pch = 20)
        }
        dev.off()
      
        pdf(paste(opt$clean_directory,"/pdfs/Crecon-", area, ".pdf", sep = ""), width = 9, height = 9)
        par(mfrow = c(5, 2))
        par(oma = c(0, 0, 0, 0))
        par(mar = c(1, 1, 1, 1))
        CreconCI <- summary(fit, pars = "Crecon", probs = c(0.025, 0.25, 0.5, 0.75, 0.975))$summary
        ind <- (Tcur - Nstep * Tstep + 1):Tcur
        CreconCI <- CreconCI[, c("2.5%", "50%", "97.5%")]
        for (i in 1:Nsample) {
          plot(t(Count[area, ]), pch = 20, ylim = c(0, max(Count[area, ind])))
          for (j in 1:3) {
            lines(ind, CreconCI[ind, j])
          }
          points(ind, Crecon_sample[area_index, ind, i], col = "red", pch = "x")
        }
        dev.off()
      }
    }
    
    days <- colnames(Count)
    rownames(Clatent_mean) <- quoted_areas
    colnames(Clatent_mean) <- days
    write.csv(Clatent_mean, paste(opt$clean_directory, "/Clatent_mean.csv", sep=""), quote = FALSE)
    rownames(Crecon_median) <- quoted_areas
    colnames(Crecon_median) <- days
    write.csv(Crecon_median, paste(opt$clean_directory, "/Crecon_median.csv", sep=""), quote = FALSE)
    for (i in 1:Nsample) {
      cc <- Clatent_sample[, , i]
      rownames(cc) <- quoted_areas
      colnames(cc) <- days
      write.csv(cc, paste(opt$clean_directory, "/Clatent_sample", i, ".csv", sep = ""), quote = FALSE)
      cc <- Crecon_sample[, , i]
      rownames(cc) <- quoted_areas
      colnames(cc) <- days
      write.csv(cc, paste(opt$clean_directory, "/Crecon_sample", i, ".csv", sep = ""), quote = FALSE)
    }

    days_modelled = days[(Tcond+1):(Tcond+Tlik)]
    Rt_percentiles = format(round(Rt_percentiles,2),nsmall=2)
    dimnames(Rt_percentiles) = list(area=quoted_areas,date=days_modelled,str_percentiles)
    df = adply(Rt_percentiles,c(2,1))
    df = df[,c(2,1,3:(num_percentiles+2))]
    write.csv(df, paste(opt$clean_directory, "/Rt.csv", sep=""), quote=FALSE, row.names=FALSE)
  })
  env
}


epiclean_cmdline_options = function(opt = epiclean_options()) {
  list(
    make_option(
      c("--task_id"),
      type="integer",
      default=1,
      help="Task ID for Slurm usage. Maps to area_index."
    ),
    make_option(
      c("--num_samples"),
      type="integer",
      default=opt$num_samples,
      help="Number of samples to output."
    ),
    make_option(
      c("--iterations"),
      type="integer",
      default=opt$iterations,
      help=paste("Number of MCMC iterations, default",opt$iterations)
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
      c("--clean_directory"), 
      type="character", 
      default=opt$clean_directory, 
      help=paste("Directory to put cleaned results in default",opt$clean_directory)
    ),

    make_option(
      c("--produce_plots"), 
      type="logical", 
      default=opt$produce_plots, 
      help=paste("Whether to produce plots; default",opt$produce_plots)
    )

  )
}

epiclean_get_cmdline_options = function(opt=epiclean_options()) {
  cmdline_opt = epiclean_cmdline_options(opt)
  opt_parser = OptionParser(option_list=cmdline_opt)
  parsed_opt = parse_args(opt_parser)
  for (o in names(parsed_opt)) {
    opt[o] = parsed_opt[o]
  }
  opt
}


