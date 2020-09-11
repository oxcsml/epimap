library(rstan)
library(geosphere)
library(optparse)

option_list = list(
  make_option(c("-s", "--spatialkernel"), type="character", default="matern12",             help="Use spatial kernel ([matern12]/matern32/matern52/exp_quad/none)"),
  make_option(c("-l", "--localkernel"),   type="character", default="local",                help="Use local kernel ([local]/none)"),
  make_option(c("-g", "--globalkernel"),  type="character", default="global",               help="Use global kernel ([global]/none)"),
  make_option(c("-m", "--metapop"),       type="character", default="radiation2,uniform,in",help="metapopulation model for inter-region cross infections (none, or comma separated list containing radiation{1,2,3},uniform,in,in_out (default is radiation2,uniform,in"),
  make_option(c("-o", "--observation"),   type="character", default="cleaned_sample",  help="observation model ([neg_binomial_{2[3]}]/poisson/cleaned_sample/cleaned_mean)"),
  make_option(c("-x", "--cleaned_sample_id"),   type="integer", default="1",  help="id of cleaned sample to use"),
  make_option(c("-c", "--chains"),        type="integer",   default=1,                      help="number of MCMC chains [4]"),
  make_option(c("-i", "--iterations"),    type="integer",   default=6000,                   help="Length of MCMC chains [6000]"),
  make_option(c("-n", "--time_steps"),    type="integer",   default=15,                      help="Number of periods to fit Rt in"),
  make_option(c("-d", "--daily_update"),  action="store_true",                              help="If True, will overide the lastest daily update of this model on compleation"),
  make_option(c("-t", "--task_id"),       type="integer",   default=0,                      help="Task ID for Slurm usage. By default, turned off [0].")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


opt$cleaned_sample_id <- cleaned_sample_id
numchains = opt$chains
numiters = opt$iterations

# If using Slurm, override other CLI options and use grid instead.
if (opt$task_id > 0) {
  grid = expand.grid(
    spatialkernel=c("matern12", "matern32", "matern52", "exp_quad", "none"), 
    metapop=c("radiation1,uniform,in", "radiation1,uniform,in_out", "radiation2,uniform,in", "radiation2,uniform,in_out", "radiation3,uniform,in", "radiation3,uniform,in_out", "uniform,in", "uniform,in_out", "none"), 
    observation=c("neg_binomial_2", "neg_binomial_3", "poisson", "cleaned_sample","cleaned_mean"),
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
#source('read_clean_data.r')
if (opt$observation == 'cleaned_sample') {
  sample_id = opt$cleaned_sample_id
  Clean_sample <- read.csv(paste('data/Clatent_sample',sample_id,'.csv',sep=''))
  print(paste('Using samples from data/Clatent_sample',sample_id,'.csv',sep=''))
} else {
  sample_id = 'mean'
  Clean_sample <- read.csv('data/Clatent_mean.csv')
  # placeholder if not using cleaned data
}

M <- opt$time_steps        # Testing with 1 time period
Tignore <- 5  # counts in most recent 7 days may not be reliable?
Tpred <- 2    # number of days held out for predictive probs eval
Tstep <- 7 # number of days to step for each time step of Rt prediction
Tlik <- M*Tstep     # number of days for likelihood to infer Rt
Tall <- Tall-Tignore  # number of days; last 7 days counts ignore; not reliable
Tcur <- Tall-Tpred       # number of days we condition on
Tcond <- Tcur-Tlik       # number of days we condition on
Tproj <- 14           # number of days to project forward

Mproj = Tproj/Tstep

AllCount <- Count
Count <- Count[,1:Tall] # get rid of ignored last days
Clean <- Clean_sample[,1:Tall] # get rid of ignored last days

days_likelihood = seq(dates[Tcond+1],by=1,length.out=Tstep*M)
days_pred_held_out = seq(dates[Tcur+1],by=1,length.out=Tpred)

print("Days used for likelihood fitting")
print(days_likelihood)
print("Days used for held out likelihood")
print(days_pred_held_out)

OBSERVATIONMODELS = list(
  'poisson' = 1,
  'neg_binomial_2' = 2,
  'neg_binomial_3' = 3,
  'cleaned_mean' = 4,
  'cleaned_sample' = 4
)
OBSERVATIONMODEL = OBSERVATIONMODELS[[opt$observation]]
if (is.null(OBSERVATIONMODEL)) {
  stop(c('Unrecognised observation option ',opt$observation));
}

# metapopulation cross-area fluxes.
METAPOPMODEL = strsplit(opt$metapop,',')[[1]]
METAPOPOPTIONS = list(
  'radiation1' = radiation_flux[,,1], # smoothed radiation model with length scale = .1 (10km)
  'radiation2' = radiation_flux[,,2], # smoothed radiation model with length scale = .2 (20km)
  'radiation3' = radiation_flux[,,3], # smoothed radiation model with length scale = .5 (50km)
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

# time steps
times = 1:M
timedist = matrix(0, M, M)
for (i in 1:M) {
  for (j in 1:M) {
    timedist[i, j] = abs(times[i] - times[j]) * Tstep
  }
}

# precompute lockdown cutoff kernel
lockdown_day = as.Date("2020-03-23")
days_lik_start = days_likelihood[seq(1, length(days_likelihood), Tstep)]
days_lik_start = vapply(days_lik_start, (function (day) as.Date(day, format="%Y-%m-%d")), double(1))
day_pre_lockdown = vapply(days_lik_start, (function (day) day < lockdown_day), logical(1))

time_corellation_cutoff = matrix(0,M,M)
for (i in 1:M) {
  for (j in 1:M) {
    time_corellation_cutoff[i, j] = !xor(day_pre_lockdown[i], day_pre_lockdown[j])
  }
}

#############################################################################################
#############################################################################################
# Main computation

Rmap_data <- list(
  N = N, 
  M = M,
  Tall = Tall,
  Tcond = Tcond,
  Tlik = Tlik,
  Tproj = Tproj,
  Tstep=Tstep,

  Count = Count,
  Clean = Clean,
  # geoloc = geoloc,
  geodist = geodist,
  timedist = timedist,
  timecorcut = time_corellation_cutoff,

  DO_METAPOP = DO_METAPOP,
  DO_IN_OUT = DO_IN_OUT,
  OBSERVATIONMODEL = OBSERVATIONMODEL,

  Tip = Tip, 
  infprofile = infprofile,
  Tdp = Tdp,
  delayprofile = delayprofile,
  F = F,
  flux = flux
)

init = list()
for (i in 1:numchains) {
  init[[i]] = list(
    gp_time_length_scale = 100.0,
    gp_space_length_scale = 2.0,
    gp_space_sigma = .01,
    global_sigma = .01,
    local_scale = .01,
    dispersion = 5.0
    )
}

runname = sprintf('Rmap-time-vary-reduce-cleaned-%s-%s-%s-%s-%s-%s-%s-%s',  
  as.character(Sys.time(),format='%Y%m%d%H%M%S'), 
  opt$spatialkernel,  
  opt$localkernel,  
  opt$globalkernel,  
  opt$metapop,  
  opt$observation, 
  opt$time_steps, 
  sample_id
)
print(runname)


# copy the stan file and put in the right kernel
stan_file_name = paste('fits/', runname, '.stan', sep='')
content = readLines(paste('stan_files/', 'Rmap-time-vary-vectorize-reduce-cleaned.stan',sep=''))
content = gsub(pattern="SPATIAL", replace=opt$spatialkernel, content)
content = gsub(pattern='TEMPORAL', replace=opt$spatialkernel, content)
content = gsub(pattern="LOCAL", replace=opt$localkernel, content)
content = gsub(pattern="GLOBAL", replace=opt$globalkernel, content)
# content = gsub(pattern="METAPOP", replace=opt$metapop, content)
#content = gsub(pattern="OBSERVATION", replace=opt$observation, content)
writeLines(content, stan_file_name)

start_time <- Sys.time()

fit <- stan(file = stan_file_name,
            data = Rmap_data, 
            iter = numiters, 
            chains = numchains,
            control = list(adapt_delta = .9))

#############################################################################################

saveRDS(fit, paste('fits/', runname, '_stanfit', '.rds', sep=''))

# fit = readRDS(paste('fits/', runname, '_stanfit', '.rds', sep=''))

end_time <- Sys.time()

print("Time to run")
print(end_time - start_time)


#############################################################################################


#################################################################
# Summary of fit
print(summary(fit, 
    pars=c("gp_space_length_scale","gp_space_sigma","gp_time_length_scale",
        "global_sigma","local_scale","dispersion",
        "Rt_all","coupling_rate","flux_probs"), 
    probs=0.5)$summary)


#################################################################
# save raw samples
save_samples = function(pars,areafirst=FALSE) {
  samples = extract(fit,pars=pars,permuted=FALSE)
  samples = samples[seq(numchains,by=numchains,to=numiters/2),,,drop=F]
  ns = numiters/2
  np = dim(samples)[3]/N
  parnames = dimnames(samples)[[3]]
  dim(samples) = c(ns,np*N)
  colnames(samples) = parnames
  if (areafirst) {
    ind = N*(rep(1:np,N)-1) + rep(1:N,each=np) 
    samples = samples[,ind,drop=F]
  } else {
    samples = samples[,,drop=F]
  }
  samples = round(samples,2)
  write.csv(samples,paste('fits/',runname,'_',pars,'_samples.csv',sep=''))
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

times = rep(c(1:M,rep(M,Mproj)), N)
places = rep(1:N, each=M+Mproj)
indicies = places + (N)*(times-1)
Rt = Rt[indicies,]

Rt <- Rt[sapply(1:(N*(M+Mproj)),function(i)rep(i,Tstep)),]
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
write.csv(df, paste('fits/', runname, '_Rt.csv', sep=''),
    row.names=FALSE,quote=FALSE)

#################################################################
# Rt exceedance probabilities
thresholds = c(.8, .9, 1.0, 1.1, 1.2, 1.5, 2.0)
numthresholds = length(thresholds)
numsamples = numchains*numiters/2
Rt <- as.matrix(fit, pars="Rt")
dim(Rt) <- c(numsamples,M,N)
Pexceedance = array(0.0,dim=c(M,N,numthresholds))
for (k in 1:M) {
  for (i in 1:N) {
    for (x in 1:numthresholds) {
      Pexceedance[k,i,x] = mean(Rt[,k,i]>thresholds[x])
    }
  }
}
Pexceedance = Pexceedance[c(1:M,rep(M,Mproj)),,]
Pexceedance <- Pexceedance[sapply(1:(M+Mproj),function(k)rep(k,Tstep)),,]
dim(Pexceedance) <- c(Tstep*(M+Mproj)*N,numthresholds)
df <- area_date_dataframe(
    quoted_areas, 
    days_all,
    provenance,
    format(round(Pexceedance,2),nsmall=2),
    c("P_08","P_09","P_10","P_11","P_12","P_15","P_20")
)
write.csv(df, paste('fits/', runname, '_Pexceed.csv', sep=''),
    row.names=FALSE,quote=FALSE)



#################################################################
# posterior predictives and projections
s <- summary(fit, pars=c("Cpred"), probs=c(0.025, .25, .5, .75, .975))$summary
Cpred <- s[,c("2.5%","25%", "50%","75%", "97.5%")]
#Cpred <- s[,c("2.5%", "50%", "97.5%")]
Cpred <- t(t(Cpred))
print(sprintf("median Cpred range: [%f, %f]",min(Cpred[,"50%"]),max(Cpred[,"50%"])))
df <- area_date_dataframe(
    quoted_areas,
    seq(dates[Tcond]+1,by=1,length.out=M*Tstep),
    rep('inferred',Tlik),
    format(round(Cpred,1),nsmall=1),
    #c("C_2_5","C_25","C_50","C_75","C_97_5")
    c("C_025","C_25","C_50","C_75","C_975")
)
write.csv(df, paste('fits/', runname, '_Cpred.csv', sep=''),
    row.names=FALSE,quote=FALSE)

# weekly counts. Includes 1 last column of actual counts among days ignored in model
Cweekly <- as.matrix(Count[,(Tcond+1):(Tcond+Tlik)])
dim(Cweekly) <- c(N,Tstep,M)
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
Cweekly <- Cweekly[sapply(1:(M+Mproj),function(k)rep(k,Tstep)),]
dim(Cweekly) <- c(N*(Tlik+Tstep*Mproj))
df <- area_date_dataframe(
    quoted_areas,
    days_all,
    provenance,
    format(round(Cweekly,1),digits=6),
    c("C_weekly")
)
write.csv(df, paste('fits/', runname, '_Cweekly.csv', sep=''),
    row.names=FALSE,quote=FALSE)



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
df <- data.frame(area = quoted_areas, logpred = logpred, provenance=rep('inferred',Tlik))
for (i in 1:Tpred)
  colnames(df)[i+1] <- sprintf('logpred_day%d',i)
write.csv(df, paste('fits/', runname, '_logpred', '.csv', sep=''),
    row.names=FALSE)

####################################################################
# pairs plot
#pdf(paste('fits/',runname,'_pairs.pdf',sep=''),width=9,height=9)
#pairs(fit, pars=c(
#    "gp_space_length_scale","gp_space_sigma","gp_time_length_scale",
#    "global_sigma","local_scale","precision")) 
#dev.off()

if (opt$daily_update) {
  dir.create('fits/latest_updates')

  runname2 = sprintf('Rmap-time-vary-reduce-cleaned-%s-%s-%s-%s-%s-%s-%s', 
    'latest',
    opt$spatialkernel, 
    opt$localkernel, 
    opt$globalkernel, 
    opt$metapop, 
    opt$observation,
    opt$time_steps,
    sample_id
  )

  file.copy(paste('fits/', runname, '_stanfit.rds', sep=''), 'fits/latest_updates', overwrite = TRUE)
  file.rename(paste('fits/latest_updates/', runname, '_stanfit.rds', sep=''), paste('fits/latest_updates/', runname2, '_stanfit.rds', sep=''))

  file.copy(paste('fits/', runname, '_Rt.csv', sep=''), 'fits/latest_updates', overwrite = TRUE)
  file.rename(paste('fits/latest_updates/', runname, '_Rt.csv', sep=''), paste('fits/latest_updates/', runname2, '_Rt.csv', sep=''))

  file.copy(paste('fits/', runname, '_Pexceed.csv', sep=''), 'fits/latest_updates', overwrite = TRUE)
  file.rename(paste('fits/latest_updates/', runname, '_Pexceed.csv', sep=''), paste('fits/latest_updates/', runname2, '_Pexceed.csv', sep=''))

  file.copy(paste('fits/', runname, '_Cpred.csv', sep=''), 'fits/latest_updates', overwrite = TRUE)
  file.rename(paste('fits/latest_updates/', runname, '_Cpred.csv', sep=''), paste('fits/latest_updates/', runname2, '_Cpred.csv', sep=''))

  file.copy(paste('fits/', runname, '_Cproj.csv', sep=''), 'fits/latest_updates', overwrite = TRUE)
  file.rename(paste('fits/latest_updates/', runname, '_Cproj.csv', sep=''), paste('fits/latest_updates/', runname2, '_Cproj.csv', sep=''))

  file.copy(paste('fits/', runname, '_Cweekly.csv', sep=''), 'fits/latest_updates', overwrite = TRUE)
  file.rename(paste('fits/latest_updates/', runname, '_Cweekly.csv', sep=''), paste('fits/latest_updates/', runname2, '_Cweekly.csv', sep=''))

  file.copy(paste('fits/', runname, '_logpred.csv', sep=''), 'fits/latest_updates', overwrite = TRUE)
  file.rename(paste('fits/latest_updates/', runname, '_logpred.csv', sep=''), paste('fits/latest_updates/', runname2, '_logpred.csv', sep=''))

}

print(runname)
