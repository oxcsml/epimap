library(rstan)
library(optparse)

option_list = list(
  make_option(c("-t", "--task_id"), type="integer", default=0, help="Task ID for Slurm usage. Maps to area_index.")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

area_index = opt$task_id
numchains = 1
numiters = 3000

options(mc.cores = min(numchains,parallel::detectCores()))
rstan_options(auto_write = TRUE)

source('read_data.r')
# counts in most recent 5 days may not be reliable
Tignore <- 0  # don't ignore for now? can ignore last 5 days of cleaned data instead?
Tall <- Tall-Tignore
Count <- Count[,1:Tall]

Nsample <- 10

Tip <- 30
infprofile <- infprofile[1:Tip]
infprofile <- infprofile/sum(infprofile)

# Gamma(5,1) days between infection and getting tested
Ttdp <- 14
Adp <- 5.0
Bdp <- 1.0
testdelayprofile <- pgamma(1:Ttdp,shape=Adp,rate=Bdp)
testdelayprofile <- testdelayprofile/testdelayprofile[Ttdp]
testdelayprofile <- testdelayprofile - c(0.0,testdelayprofile[1:(Ttdp-1)])

# Case only reported a few days after testing, 
# no result delay truncation
Trdp <- 5
resultdelaydecay = .5
resultdelaystrength = 5
# Geometric(.5) delay distribution.
# Trdp <- 5
# resultdelayprofile <- .5^seq(1,by=1,length.out=Trdp)
# resultdelayprofile <- resultdelayprofile / sum(resultdelayprofile)


Tstep <- 7
#Nstep <- floor((Tall-max(Tip,Ttdp)) / Tstep)
Nstep <- 18 # about 4 months
Tlik <- Nstep*Tstep
Tcond <- Tall-Tlik


area = areas[area_index]
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
  Ttdp = Ttdp,
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

fit <- stan(file = 'stan_files/Rmap-clean.stan',
            data = Rmap_clean_data, 
            init = init,
            iter = numiters, 
            chains = numchains,
            control = list(adapt_delta = .9))

end_time <- Sys.time()

print("Time to run")
print(end_time - start_time)

saveRDS(fit, paste('local-cleaned/stanfit-',area,'.rds',sep=''))

#################################################################
# Summary of fit
print(
    summary(fit, 
        pars=c("mu","sigma","alpha","phi_latent","phi_observed","xi","Noutliers","meandelay","resultdelayprofile","Rt"),
        probs=c(0.5)
    )$summary
)
