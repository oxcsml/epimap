library(rstan)
library(geosphere)
library(optparse)

option_list = list(
  make_option(c("-k", "--kernel"), type="character", default="matern12", 
              help="kernel to use in the spatial prior GP ([matern12]/matern32/matern52/exp_quad/none)"),
  make_option(c("-m", "--metapop"), type="character", default="none", 
              help="metapopulation model for inter-region cross infections ([uniform1]/uniform2/none)"),
  make_option(c("-l", "--likelihood"), type="character", default="negative_binomial", 
              help="likelihood model ([negative_binomial]/poisson)"),
  make_option(c("-c", "--chains"), type="integer", default=4,
              help="number of MCMC chains [4]"),
  make_option(c("-i", "--iterations"), type="integer", default=4000,
              help="Length of MCMC chains [4000]")
); 


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

numchains = opt$chains
numiters = opt$iterations

options(mc.cores = min(numchains,parallel::detectCores()))
rstan_options(auto_write = TRUE)

infprofile <- read.csv("data/serial_interval.csv")$fit

uk_cases <- read.csv("data/uk_cases.csv")

metadata <- read.csv("data/metadata.csv")

N <- 149                 # number of regions in England only
D <- 100                 # infection profile number of days
Tignore <- 7             # counts in most recent 7 days may not be reliable?
Tpred <- 7               # number of days held out for predictive probs eval
Tlik <- 7                # number of days for likelihood to infer Rt
Tall <- ncol(uk_cases)-2-Tignore  # number of days; last 7 days counts ignore; not reliable
Tcond <- Tall-Tlik-Tpred       # number of days we condition on
Tproj <- 21              # number of days to project forward


Count <- uk_cases[1:N,3:(Tall+2)]

geoloc <- matrix(0, N, 2)
geodist <- matrix(0, N, N)

region_names <- metadata$AREA
longitudes <- metadata$LONG
latitudes <- metadata$LAT

for (i in 1:N) {
  region_name <- uk_cases[i,2]
  j <- grep(sprintf('^%s$',region_name), region_names)
  if (length(j) >= 1) {                 # just use first match!!
    if (length(j) > 1) 
      print(sprintf("Found regions %s, using first", paste(region_names[j]), collapse=","))
    geoloc[i, 1] = longitudes[j[1]]
    geoloc[i, 2] = latitudes[j[1]]
  } else {
    print(sprintf("Cannot find region '%s'",region_name))
    for (r in 1:length(region_names)) {
      if (length(grep(region_names[r], region_name))>0) {
        geoloc[i, 1] = longitudes[r]
        geoloc[i, 2] = latitudes[r]
        print(sprintf("...found region '%s'",region_names[r]))
      }
    }
  }
}
  
for (i in 1:N) {
  for (j in i:N) {
    # distance between two points on an ellipsoid (default is WGS84 ellipsoid), in units of 100km
    geodist[i, j] = distGeo(geoloc[i, 1:2], geoloc[j, 1:2]) / 100000
    geodist[j, i] = geodist[i, j]
  }
}


Rmap_data <- list(N = N, D = D, 
                 Tall = Tall, 
                 Tcond = Tcond, 
                 Tlik = Tlik, 
                 Tproj = Tproj, 
                 Count = Count,  
                 geoloc = geoloc,
                 geodist = geodist,
                 infprofile = infprofile)

# copy the stan file and put in the right kernel
runname = sprintf('Rmap-%s-%s-%s', opt$kernel, opt$metapop, opt$likelihood)
stan_file_name = paste('fits/', runname, '.stan', sep='')
stan_file_content = readLines(paste('stan_files/', 'Rmap.stan',sep=''))
stan_file_content = gsub(pattern="KERNEL", replace=opt$kernel, x=
                    gsub(pattern="METAPOP", replace=opt$metapop, x=
                    gsub(pattern="LIKELIHOOD", replace=opt$likelihood, x=
                    stan_file_content)))
writeLines(stan_file_content, stan_file_name)

fit <- stan(file = stan_file_name,
            data = Rmap_data, 
            iter = numiters, 
            chains = numchains,
            control = list(adapt_delta = .9))
# print(fit)

print(summary(fit, 
    pars=c("Ravg","gp_length_scale","gp_sigma","global_sigma","local_sigma","dispersion","coupling_rate"), 
    probs=c(0.025, 0.5, 0.975))$summary)


s <- summary(fit, pars="Rt", probs=c(0.025, .5, .975))$summary
Rt <- s[,c("2.5%","50%","97.5%")]
Rt <- t(t(Rt))

print(sprintf("median Rt range: [%f, %f]",min(Rt[,2]),max(Rt[,2])))

s <- summary(fit, pars="Ppred", probs=c(0.025, .5, .975))$summary
Ppred <- s[,"mean"]
logpred <- log(Ppred)
dim(logpred) <- c(Tpred,N)
logpred <- t(logpred)
print(sprintf("mean log predictives = %f",mean(logpred)))


s <- summary(fit, pars="Cproj", probs=c(0.025, .5, .975))$summary
Cproj <- s[,c("2.5%","50%","97.5%")]
Cproj <- t(t(Cproj))
Cprojlower <- Cproj[,1]
Cprojmedian <- Cproj[,2]
Cprojupper <- Cproj[,3]
dim(Cprojlower) <- c(Tproj,N)
dim(Cprojmedian) <- c(Tproj,N)
dim(Cprojupper) <- c(Tproj,N)
Cprojlower <- t(Cprojlower)
Cprojmedian <- t(Cprojmedian)
Cprojupper <- t(Cprojupper)

print(sprintf("median Cproj range: [%f, %f]",min(Cproj[,2]),max(Cproj[,2])))

df <- data.frame(area = uk_cases[1:N,2], logpred = logpred)
colnames(df)[1] <- "area"
for (i in 1:Tpred)
  colnames(df)[i+1] <- sprintf('logpred_day%d',i)
write.csv(df, paste('fits/', 'logpred_', runname, '.csv', sep=''),row.names=FALSE)

df <- data.frame(area = uk_cases[1:N,2], Rt = Rt, Cproj = Cproj)
colnames(df) <- c("area","Rtlower","Rtmedian","Rtupper","Cprojlower","Cprojmedian","Cprojupper")
write.csv(df, paste('fits/', 'RtCproj_', runname, '.csv', sep=''),row.names=FALSE)

saveRDS(fit, paste('fits/', 'stanfit_', runname, '.rds', sep=''))

print(runname)

