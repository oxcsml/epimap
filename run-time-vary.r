library(rstan)
library(geosphere)
library(optparse)

option_list = list(
  make_option(c("-s", "--spatialkernel"), type="character", default="matern12",             help="Use spatial kernel ([matern12]/matern32/matern52/exp_quad/none)"),
  make_option(c("-l", "--localkernel"),   type="character", default="local",                help="Use local kernel ([local]/none)"),
  make_option(c("-g", "--globalkernel"),  type="character", default="global",               help="Use global kernel ([global]/none)"),
  make_option(c("-m", "--metapop"),       type="character", default="radiation2_uniform_in",help="metapopulation model for inter-region cross infections (uniform_in{_out}/[radiation{1[2]3}_uniform_in{_out}]/none)"),
  make_option(c("-o", "--observation"),   type="character", default="negative_binomial_3",  help="observation model ([negative_binomial_{2[3]}]/poisson)"),
  make_option(c("-c", "--chains"),        type="integer",   default=4,                      help="number of MCMC chains [4]"),
  make_option(c("-i", "--iterations"),    type="integer",   default=6000,                   help="Length of MCMC chains [6000]"),
  make_option(c("-n", "--time_steps"),    type="integer",   default=6,                      help="Number of periods to fit Rt in"),
  make_option(c("-t", "--task_id"),       type="integer",   default=0,                      help="Task ID for Slurm usage. By default, turned off [0].")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

numchains = opt$chains
numiters = opt$iterations

# If using Slurm, override other CLI options and use grid instead.
if (opt$task_id > 0) {
  grid = expand.grid(
    spatialkernel=c("matern12", "matern32", "matern52", "exp_quad", "none"), 
    metapop=c("radiation1_uniform_in", "radiation1_uniform_in_out", "radiation2_uniform_in", "radiation2_uniform_in_out", "radiation3_uniform_in", "radiation3_uniform_in_out", "uniform_in", "uniform_in_out", "none"), 
    observation=c("negative_binomial_2", "negative_binomial_3", "poisson"),
    localkernel=c("local","none"),
    globalkernel=c("global","none")
  )
  grid = sapply(grid, as.character)
  update = as.list(grid[opt$task_id, ]) 
  for (name in names(update)){
    opt[name] = update[name]
  }
}

options(mc.cores = min(numchains,parallel::detectCores()))
rstan_options(auto_write = TRUE)

source('read_data.r')
source('read_radiation_fluxes.r')

M <- opt$time_steps        # Testing with 1 time period
Tignore <- 4  # counts in most recent 7 days may not be reliable?
Tpred <- 3    # number of days held out for predictive probs eval
Tlik <- 7     # number of days for likelihood to infer Rt
Tstep <- Tlik # number of days to step for each time step of Rt prediction
Tall <- Tall-Tignore  # number of days; last 7 days counts ignore; not reliable
Tcond <- Tall-Tlik-Tpred       # number of days we condition on
Tproj <- 7            # number of days to project forward

Count <- Count[,1:Tall] # get rid of ignored last days

days_likelihood = seq(dates[Tcond - (M-1)*Tlik +1],by=1,length.out=Tlik*M)
days_pred_held_out = seq(dates[Tcond+Tlik+1],by=1,length.out=Tpred)

print("Days used for likelihood fitting")
print(days_likelihood)
print("Days used for held out likelihood")
print(days_pred_held_out)

# metapopulation cross-area fluxes.
if (opt$metapop == 'radiation1_uniform_in' || 
    opt$metapop == 'radiation1_uniform_in_out') {
  do_metapop = 1
  if (opt$metapop == 'radiation1_uniform_in' ) {
    do_in_out = 0
  } else {
    do_in_out = 1
  }
  flux = list()
  flux[[1]] = radiation_flux[,,1] # ls=.1
  flux[[2]] = matrix(1.0/N,N,N) # uniform cross-area infections
  F = length(flux)
} else if (opt$metapop == 'radiation2_uniform_in' || 
           opt$metapop == 'radiation2_uniform_in_out') {
  do_metapop = 1
  if (opt$metapop == 'radiation2_uniform_in' ) {
    do_in_out = 0
  } else {
    do_in_out = 1
  }
  flux = list()
  flux[[1]] = radiation_flux[,,2] # ls=.1
  flux[[2]] = matrix(1.0/N,N,N) # uniform cross-area infections
  F = length(flux)
} else if (opt$metapop == 'radiation3_uniform_in' || 
           opt$metapop == 'radiation3_uniform_in_out') {
  do_metapop = 1
  if (opt$metapop == 'radiation3_uniform_in' ) {
    do_in_out = 0
  } else {
    do_in_out = 1
  }
  flux = list()
  flux[[1]] = radiation_flux[,,3] # ls=.1
  flux[[2]] = matrix(1.0/N,N,N) # uniform cross-area infections
  F = length(flux)
} else if (opt$metapop == 'uniform_in' || 
           opt$metapop == 'uniform_in_out') {
  do_metapop = 1
  if (opt$metapop == 'uniform_in' ) {
    do_in_out = 0
  } else {
    do_in_out = 1
  }
  flux = list()
  flux[[1]] = matrix(1.0/N,N,N) # uniform cross-area infections
  F = length(flux)
} else if (opt$metapop == 'none') {
  do_metapop = 0
  do_in_out = 0
  flux = array(0,dim=c(0,N,N));
  F = 0
} else {
  stop(c('Unrecognised metapop option ',opt$metapop));
}

times = (1:M) * Tstep
timedist = matrix(0, M, M)
for (i in 1:M) {
  for (j in 1:M) {
    timedist[i, j] = abs(times[i] - times[j]) * Tstep
  }
}

# precompute lockdown cutoff kernel
lockdown_day = as.Date("2020-03-23")
days_lik_start = days_likelihood[seq(1, length(days_likelihood), Tlik)]
days_lik_start = vapply(days_lik_start, (function (day) as.Date(day, format="%Y-%m-%d")), double(1))
day_pre_lockdown = vapply(days_lik_start, (function (day) day < lockdown_day), logical(1))

time_corellation_cutoff = matrix(0,M,M)
for (i in 1:M) {
  for (j in 1:M) {
    time_corellation_cutoff[i, j] = !xor(day_pre_lockdown[i], day_pre_lockdown[j])
  }
}

Rmap_data <- list(
  N = N, 
  M = M,
  D = D, 
  Tall = Tall,
  Tcond = Tcond,
  Tlik = Tlik,
  Tproj = Tproj,
  Tstep=Tstep,
  Count = Count,
  # geoloc = geoloc,
  geodist = geodist,
  timedist = timedist,
  timecorcut = time_corellation_cutoff,
  do_metapop = do_metapop,
  do_in_out = do_in_out,
  F = F,
  flux = flux,
  infprofile = infprofile
  # local_sd = opt$local_sd,
  # global_sd = opt$global_sd,
  # gp_sd = opt$gp_sd,
  # gp_length_scale_sd = opt$gp_length_scale_sd
)

runname = sprintf('Rmap-time-vary-%s-%s-%s-%s-%s-%s-%s', 
  as.character(Sys.time(),format='%Y%m%d%H%M%S'),
  opt$spatialkernel, 
  opt$localkernel, 
  opt$globalkernel, 
  opt$metapop, 
  opt$observation,
  opt$time_steps
)
print(runname)


# copy the stan file and put in the right kernel
stan_file_name = paste('fits/', runname, '.stan', sep='')
content = readLines(paste('stan_files/', 'Rmap-time-vary-vectorize.stan',sep=''))
content = gsub(pattern="SPATIAL", replace=opt$spatialkernel, content)
content = gsub(pattern='TEMPORAL', replace=opt$spatialkernel, content)
content = gsub(pattern="LOCAL", replace=opt$localkernel, content)
content = gsub(pattern="GLOBAL", replace=opt$globalkernel, content)
# content = gsub(pattern="METAPOP", replace=opt$metapop, content)
content = gsub(pattern="OBSERVATION", replace=opt$observation, content)
writeLines(content, stan_file_name)

start_time <- Sys.time()

fit <- stan(file = stan_file_name,
            data = Rmap_data, 
            iter = numiters, 
            chains = numchains,
            control = list(adapt_delta = .9))
saveRDS(fit, paste('fits/', runname, '_stanfit', '.rds', sep=''))

end_time <- Sys.time()

print("Time to run")
print(end_time - start_time)


# fit = readRDS(paste('fits/', runname, '_stanfit', '.rds', sep=''))


#################################################################
# Summary of fit
print(summary(fit, 
    pars=c("R0","gp_space_length_scale","gp_space_sigma","gp_time_length_scale","global_sigma","local_scale","precision","coupling_rate"), 
    probs=c(0.025, 0.25, 0.5, 0.75, 0.975))$summary)
  
#################################################################
area_date_dataframe <- function(areas,dates,data,data_names) {
  numareas <- length(areas)
  numdates <- length(dates)
  dates <- rep(dates,numareas)
  dim(dates) <- c(numareas*numdates)
  areas <- rep(areas,numdates)
  dim(areas) <- c(numareas,numdates)
  areas <- t(areas)
  dim(areas) <- c(numareas*numdates)
  df <- data.frame(area=areas,Date=dates,data=data)
  colnames(df)[3:ncol(df)] <- data_names
  df
}

#################################################################
# Rt posterior
#s <- summary(fit, pars="Rt", probs=c(0.025, .1, .2, 0.25, .3, .4, .5, .6, .7, 0.75, .8, .9, .975))$summary
#Rt <- s[,c("2.5%","10%", "20%", "25%", "30%", "40%", "50%", "60%", "70%", "75%", "80%", "90%","97.5%")]
s <- summary(fit, pars="Rt", probs=c(.1, .2, .3, .4, .5, .6, .7, .8, .9))$summary
Rt <- s[,c("10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%")]

times = rep(1:M, N)
places = rep(1:N, each=M)
indicies = places + (N)*(times-1)
Rt = Rt[indicies,]

Rt <- Rt[sapply(1:(N*M),function(i)rep(i,Tlik)),]
print(sprintf("median Rt range: [%f, %f]",min(Rt[,"50%"]),max(Rt[,"50%"])))
df <- area_date_dataframe(
    quoted_areas, 
    #dates[Tcond+1],
    days_likelihood,
    format(Rt,digits=2),
    c("Rt_10","Rt_20","Rt_30","Rt_40","Rt_50","Rt_60","Rt_70","Rt_80","Rt_90")
    #c("Rt_2_5","Rt_10","Rt_20","Rt_25","Rt_30","Rt_40","Rt_50",
    #  "Rt_60","Rt_70","Rt_75","Rt_80","Rt_90","Rt_97_5")
)
df <- df[,c(1,3,4,5,6,7,8,9,10,11,2)]
write.csv(df, paste('fits/', runname, '_Rt.csv', sep=''),
    row.names=FALSE,quote=FALSE)

#################################################################
# projections
s <- summary(fit, pars="Cproj", probs=c(0.025, .25, .5, .75, .975))$summary
#Cproj <- s[,c("2.5%","25%", "50%","75%", "97.5%")]
Cproj <- s[,c("2.5%", "50%", "97.5%")]
Cproj <- t(t(Cproj))
print(sprintf("median Cproj range: [%f, %f]",min(Cproj[,"50%"]),max(Cproj[,"50%"])))
df <- area_date_dataframe(
    quoted_areas,
    seq(dates[Tcond+Tlik]+1,by=1,length.out=Tproj),
    format(Cproj,digits=2),
    #c("C_2_5","C_25","C_50","C_75","C_97_5")
    c("C_lower","C_median","C_upper")
)
write.csv(df, paste('fits/', runname, '_Cproj.csv', sep=''),
    row.names=FALSE,quote=FALSE)

#################################################################
# predictive probabilities
s <- summary(fit, pars="Ppred", probs=c(0.025, .5, .975))$summary
Ppred <- s[,"mean"]
logpred <- log(Ppred)
dim(logpred) <- c(Tpred,N)
logpred <- t(logpred)
print(sprintf("mean log predictives = %f",mean(logpred)))
df <- data.frame(area = quoted_areas, logpred = logpred)
for (i in 1:Tpred)
  colnames(df)[i+1] <- sprintf('logpred_day%d',i)
write.csv(df, paste('fits/', runname, '_logpred', '.csv', sep=''),
    row.names=FALSE)

####################################################################
# pairs plot
pdf(paste('fits/',runname,'_pairs.pdf',sep=''),width=9,height=9)
pairs(fit, pars=c(
    "R0","gp_space_length_scale","gp_space_sigma","gp_time_length_scale",
    "global_sigma","local_scale","precision","coupling_rate")) 
dev.off()

print(runname)
