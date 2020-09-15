library(rstan)

source("read_data.r")

numiters <- 3000 # NOTE: this should match the value used in clean_areas.r. TODO: Should we put this variable in a common file from which to import?

# counts in most recent 5 days may not be reliable
Tignore <- 5 # don't ignore for now? can ignore last 5 days of cleaned data instead?
Tall <- Tall - Tignore
Count <- Count[, 1:Tall]

Nsample <- 10

Tip <- 30
infprofile <- infprofile[1:Tip]
infprofile <- infprofile / sum(infprofile)

# Gamma(5,1) days between infection and getting tested
Ttdp <- 14
Adp <- 5.0
Bdp <- 1.0
testdelayprofile <- pgamma(1:Ttdp, shape = Adp, rate = Bdp)
testdelayprofile <- testdelayprofile / testdelayprofile[Ttdp]
testdelayprofile <- testdelayprofile - c(0.0, testdelayprofile[1:(Ttdp - 1)])

# Case only reported a few days after testing,
# no result delay truncation
Trdp <- 1
resultdelayprofile <- array(1)
# Geometric(.5) delay distribution.
# Trdp <- 5
# resultdelayprofile <- .5^seq(1,by=1,length.out=Trdp)
# resultdelayprofile <- resultdelayprofile / sum(resultdelayprofile)


Tstep <- 7
# Nstep <- floor((Tall-max(Tip,Ttdp)) / Tstep)
Nstep <- 18 # about 4 months
Tlik <- Nstep * Tstep
Tcond <- Tall - Tlik


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

  fit <- readRDS(paste("local-cleaned/stanfit-", area, ".rds", sep = ""))

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
  pdf(paste("local-cleaned/pairs-", area, ".pdf", sep = ""), width = 9, height = 9)
  pairs(fit, pars = c("mu", "sigma", "alpha", "phi_latent", "phi_observed"))
  dev.off()

  pdf(paste("local-cleaned/Clatent-", area, ".pdf", sep = ""), width = 9, height = 9)
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

  pdf(paste("local-cleaned/Crecon-", area, ".pdf", sep = ""), width = 9, height = 9)
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
write.csv(Clatent_mean, "data/Clatent_mean.csv", quote = FALSE)
rownames(Crecon_median) <- quoted_areas
colnames(Crecon_median) <- days
write.csv(Crecon_median, "data/Crecon_median.csv", quote = FALSE)
for (i in 1:Nsample) {
  cc <- Clatent_sample[, , i]
  rownames(cc) <- quoted_areas
  colnames(cc) <- days
  write.csv(cc, paste("data/Clatent_sample", i, ".csv", sep = ""), quote = FALSE)
  cc <- Crecon_sample[, , i]
  rownames(cc) <- quoted_areas
  colnames(cc) <- days
  write.csv(cc, paste("data/Crecon_sample", i, ".csv", sep = ""), quote = FALSE)
}