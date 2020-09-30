library(rstan)
library(optparse)
source("dataprocessing/read_data.r")

epiclean_options = function(
  num_samples     = 10,
  iterations      = 3000,
  num_steps       = NULL,
  days_per_step   = 7,
  days_ignored    = 0,
  data_directory  = "data/",
  clean_directory = "results/default"
) {
  as.list(environment())
}

epiclean = function(area_index = 0, opt = epiclean_options()) {
  env = new.env(parent=globalenv())
  env$area_index = area_index
  env$opt = opt
  Rmap_read_data(env)
  with(env,{

    numchains = 1
    numiters = opt$iterations

    options(mc.cores = min(numchains,parallel::detectCores()))

    # counts in most recent 5 days may not be reliable
    Tignore <- opt$days_ignored  # don't ignore for now? can ignore last 5 days of cleaned data instead?
    Tstep <- opt$days_per_step
    Tall <- Tall-Tignore
    if (is.null(opt$num_steps)) {
      Nstep <- floor((Tall-length(testdelayprofile))/Tstep)
      print(paste("Nstep = ",Nstep))
    } else {
      Nstep <- opt$num_steps
    }
    Tlik <- Nstep*Tstep
    Tcond <- Tall-Tlik

    Nsample <- opt$num_samples

    # Case only reported a few days after testing, 
    # no result delay truncation
    Trdp <- 5
    resultdelaydecay = .5
    resultdelaystrength = 5

    area = areas[area_index]
    Count <- AllCount[area,1:Tall]
    print(area)

    # ---------------------------------------------------------------------------- #
    #                               Main computation                               #
    # ---------------------------------------------------------------------------- #

    Rmap_clean_data <- list(
      Tall = Tall,
      Tstep = Tstep, 
      Nstep = Nstep,
      Count = Count[area,],
      Tip = Tip,
      infprofile = infprofile,
      Ttdp = length(testdelayprofile),
      testdelayprofile = testdelayprofile,
      Trdp = Trdp,
      resultdelaydecay = resultdelaydecay,
      resultdelaystrength = resultdelaystrength,
      mu_scale = 0.5,
      sigma_scale = 0.5,
      alpha_scale = 0.1,
      phi_latent_scale = 1.0,
      phi_observed_scale = 5.0,
      outlier_prob_threshold = .95,
      outlier_count_threshold = 2,
      xi_scale = 0.1,
      reconstruct_infections = TRUE
    )
    
    init = list()
    init[[1]] = list(
      mu = 0.0,
      sigma = 0.01,
      alpha1 = .95,
      phi_latent = 1.0,
      phi_observed = 5.0,
      xi = .01,
      'Reta[1]' = 0.0,
      'Reta[2]' = 0.0,
      'Reta[3]' = 0.0,
      'Reta[4]' = 0.0,
      'Reta[5]' = 0.0,
      'Reta[6]' = 0.0,
      'Reta[7]' = 0.0,
      'Reta[8]' = 0.0,
      'Reta[9]' = 0.0,
      'Reta[10]' = 0.0,
      'Reta[11]' = 0.0,
      'Reta[12]' = 0.0,
      'Reta[13]' = 0.0,
      'Reta[14]' = 0.0,
      'Reta[15]' = 0.0
    )
    
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
    
    saveRDS(fit, paste(opt$clean_directory, '/stanfits/',area,'.rds',sep=''))
    
    #################################################################
    # Summary of fit
    print(summary(fit, 
      pars=c("mu","sigma","alpha","phi_latent","phi_observed","xi","Noutliers","meandelay","resultdelayprofile","Rt"),
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

    numiters <- opt$iterations 

    Tignore <- opt$days_ignored  # don't ignore for now? can ignore last 5 days of cleaned data instead?
    Tstep <- opt$days_per_step
    Tall <- Tall-Tignore
    if (is.null(opt$num_steps)) {
      Nstep <- floor((Tall-length(testdelayprofile))/Tstep)
      print(paste("Nstep = ",Nstep))
    } else {
      Nstep <- opt$num_steps
    }
    Tlik <- Nstep*Tstep
    Tcond <- Tall-Tlik

    Nsample <- opt$num_samples

    Count <- AllCount[, 1:Tall]


    # Initialize arrays
    Clatent_sample <- array(0, c(N, Tall, Nsample))
    Clatent_mean <- array(0, c(N, Tall))
    Crecon_sample <- array(0, c(N, Tall, Nsample))
    Crecon_median <- array(0, c(N, Tall))
    Clatent_mean[, 1:Tcond] <- as.matrix(Count[, 1:Tcond])
    for (i in 1:Nsample) {
      Clatent_sample[, 1:Tcond, i] <- as.matrix(Count[, 1:Tcond])
      Crecon_sample[, 1:Tcond, i] <- as.matrix(Count[, 1:Tcond])
    }
    
    # Loop over areas, loading area RDS files and filling the arrays
    for (area_index in 1:N) {
      area <- areas[area_index]
      print(area)
    
      fit <- readRDS(paste(opt$clean_directory, '/stanfits/',area,'.rds',sep=''))
    
      skip <- numiters / 2 / Nsample
      ####################################################################
      Clatent_s <- extract(fit, pars = "Clatent", permuted = FALSE)
      Clatent_s <- Clatent_s[seq(skip, by = skip, length.out = Nsample), , ]
      dim(Clatent_s) <- c(Nsample, Tall)
      Clatent_s <- t(Clatent_s)
      dim(Clatent_s) <- c(1, Tall, Nsample)
      Clatent_sample[area_index, , ] <- Clatent_s
      Clatent_m <- summary(fit, pars = "Clatent", probs = c(0.5))$summary
      Clatent_m <- t(as.matrix(Clatent_m[, "mean"]))
      Clatent_mean[area_index, ] <- Clatent_m

      ####################################################################
      Crecon_s <- extract(fit, pars = "Crecon", permuted = FALSE)
      Crecon_s <- Crecon_s[seq(skip, by = skip, length.out = Nsample), , ]
      Crecon_s <- t(Crecon_s)
      dim(Crecon_s) <- c(1, Tall, Nsample)
      Crecon_sample[area_index, , ] <- Crecon_s
      Crecon_m <- summary(fit, pars = "Crecon", probs = c(0.5))$summary
      Crecon_m <- t(as.matrix(Crecon_m[, "50%"]))
      Crecon_median[area_index, ] <- round(Crecon_m)


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
      ind <- (Tcond + 1):Tall
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
      ind <- (Tall - Nstep * Tstep + 1):Tall
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
  })
  env
}
