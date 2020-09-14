
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
Clean_sample <- read.csv('data/Clatent_mean.csv')

M <- opt$time_steps        # Testing with 1 time period
Tignore <- 6  # counts in most recent 7 days may not be reliable?
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


runs_latent = c(
  'Rmap-time-vary-reduce-cleaned-20200910225318-matern12-local-global-radiation2,uniform,in-cleaned_sample-15',
  'Rmap-time-vary-reduce-cleaned-20200910225535-matern12-local-global-radiation2,uniform,in-cleaned_sample-15',
  'Rmap-time-vary-reduce-cleaned-20200910230223-matern12-local-global-radiation2,uniform,in-cleaned_sample-15',
  'Rmap-time-vary-reduce-cleaned-20200910224547-matern12-local-global-radiation2,uniform,in-cleaned_sample-15',
  'Rmap-time-vary-reduce-cleaned-20200910231253-matern12-local-global-radiation2,uniform,in-cleaned_sample-15',
  'Rmap-time-vary-reduce-cleaned-20200910211954-matern12-local-global-radiation2,uniform,in-cleaned_sample-15',
  'Rmap-time-vary-reduce-cleaned-20200910212148-matern12-local-global-radiation2,uniform,in-cleaned_sample-15',
  'Rmap-time-vary-reduce-cleaned-20200910205938-matern12-local-global-radiation2,uniform,in-cleaned_sample-15',
  'Rmap-time-vary-reduce-cleaned-20200910212111-matern12-local-global-radiation2,uniform,in-cleaned_sample-15'
)
runs_recon = c(
  'Rmap-time-vary-reduce-cleaned-matern12-local-global-radiation2,uniform,in-cleaned_recon_sample_1-steps_15',
  'Rmap-time-vary-reduce-cleaned-matern12-local-global-radiation2,uniform,in-cleaned_recon_sample_2-steps_15',
  'Rmap-time-vary-reduce-cleaned-matern12-local-global-radiation2,uniform,in-cleaned_recon_sample_3-steps_15',
  'Rmap-time-vary-reduce-cleaned-matern12-local-global-radiation2,uniform,in-cleaned_recon_sample_4-steps_15',
  'Rmap-time-vary-reduce-cleaned-matern12-local-global-radiation2,uniform,in-cleaned_recon_sample_5-steps_15',
  'Rmap-time-vary-reduce-cleaned-matern12-local-global-radiation2,uniform,in-cleaned_recon_sample_6-steps_15',
  'Rmap-time-vary-reduce-cleaned-matern12-local-global-radiation2,uniform,in-cleaned_recon_sample_7-steps_15',
  'Rmap-time-vary-reduce-cleaned-matern12-local-global-radiation2,uniform,in-cleaned_recon_sample_8-steps_15',
  'Rmap-time-vary-reduce-cleaned-matern12-local-global-radiation2,uniform,in-cleaned_recon_sample_9-steps_15',
  'Rmap-time-vary-reduce-cleaned-matern12-local-global-radiation2,uniform,in-cleaned_recon_sample_10-steps_15'
)
runs = runs_recon
numruns = length(runs)

allrunname = 'Rmap-time-vary-reduce-cleaned-allsamples-matern12-local-global-radiation2,uniform,in-cleaned_sample-15'

load_samples = function(pars) {
  samples = do.call(rbind, lapply(1:length(runs), function(i) {
    read.csv(paste('fits/',runs[i],'_',pars,'_samples.csv',sep=''))
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

Rt = Rt[sapply(1:N,function(i)rep((i-1)*M+c(1:M,rep(M,Mproj)),each=Tstep)),]
df <- area_date_dataframe(
    quoted_areas, 
    days_all,
    provenance,
    format(round(Rt,2),nsmall=2),
    #c("Rt_10","Rt_20","Rt_30","Rt_40","Rt_50","Rt_60","Rt_70","Rt_80","Rt_90")
    c("Rt_2_5","Rt_10","Rt_20","Rt_25","Rt_30","Rt_40","Rt_50",
      "Rt_60","Rt_70","Rt_75","Rt_80","Rt_90","Rt_97_5")
)
write.csv(df, paste('fits/', allrunname, '_Rt.csv', sep=''),
    row.names=FALSE,quote=FALSE)

#################################################################
# Rt exceedance probabilities
thresholds = c(.8, .9, 1.0, 1.1, 1.2, 1.5, 2.0)
numthresholds = length(thresholds)
numsamples = numruns * numiters/2
Rt <- as.matrix(Rt_samples)
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
write.csv(df, paste('fits/', allrunname, '_Pexceed.csv', sep=''),
    row.names=FALSE,quote=FALSE)



#################################################################
# posterior predictives and projections
Cpred = t(apply(Cpred_samples,2,quantile,
    probs=c(0.025, 0.25, .5, 0.75, .975)
))

df <- area_date_dataframe(
    quoted_areas,
    seq(dates[Tcond]+1,by=1,length.out=M*Tstep),
    rep('inferred',Tlik),
    format(round(Cpred,1),nsmall=1),
    #c("C_2_5","C_25","C_50","C_75","C_97_5")
    c("C_025","C_25","C_50","C_75","C_975")
)
write.csv(df, paste('fits/', allrunname, '_Cpred.csv', sep=''),
    row.names=FALSE,quote=FALSE)

####################################################################################
# weekly counts. Includes 1 last column of actual counts among days ignored in model
Cweekly <- as.matrix(Count[,(Tcond+1):(Tcond+Tlik)])
dim(Cweekly) <- c(N,Tstep,M)
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
Cweekly <- Cweekly[sapply(1:(M+Mproj),function(k)rep(k,Tstep)),]
dim(Cweekly) <- c(N*(Tlik+Tstep*Mproj))
df <- area_date_dataframe(
    quoted_areas,
    days_all,
    provenance,
    format(Cweekly,digits=3),
    c("C_weekly")
)
write.csv(df, paste('fits/', allrunname, '_Cweekly.csv', sep=''),
    row.names=FALSE,quote=FALSE)



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
write.csv(df, paste('fits/', allrunname, '_Cproj.csv', sep=''),
    row.names=FALSE,quote=FALSE)


