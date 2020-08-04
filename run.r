library(rstan)
library(geosphere)
library(optparse)

option_list = list(
  make_option(c("-s", "--spatialkernel"), type="character",default="matern12",   help="Use spatial kernel ([matern12]/matern32/matern52/exp_quad/none)"),
  make_option(c("-l", "--localkernel"),   type="character",default="local",    help="Use local kernel ([local]/none)"),
  make_option(c("-g", "--globalkernel"),  type="character",default="global",    help="Use global kernel ([global]/none)"),
  make_option(c("-m", "--metapop"),       type="character",default="radiation_uniform_in",   help="metapopulation model for inter-region cross infections (uniform_out/uniform_in/radiation_out/radiation_in/[radiation_uniform_in]/none)"),
  make_option(c("-o", "--observation"),   type="character",default="negative_binomial_3", help="observation model ([negative_binomial]/poisson)"),
  make_option(c("-c", "--chains"),        type="integer",  default=4,        help="number of MCMC chains [4]"),
  make_option(c("-i", "--iterations"),    type="integer",  default=6000,     help="Length of MCMC chains [4000]"),
  make_option(c("-t", "--task_id"), type="integer", default=0,               help="Task ID for Slurm usage. By default, turned off [0].")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

numchains = opt$chains
numiters = opt$iterations

# If using Slurm, override other CLI options and use grid instead.
if (opt$task_id > 0) {
  grid = expand.grid(
    spatialkernel=c("matern12", "matern32", "matern52", "exp_quad", "none"), 
    localkernel=c("local","none"),
    globalkernel=c("global","none"),
    metapop=c("uniform1", "uniform2", "none"), 
    observation=c("negative_binomial", "poisson")
  )
  grid = sapply(grid, as.character)
  opt = as.list(grid[opt$task_id, ])  
}

options(mc.cores = min(numchains,parallel::detectCores()))
rstan_options(auto_write = TRUE)

infprofile <- read.csv("data/serial_interval.csv")$fit

uk_cases <- read.csv("data/uk_cases.csv")
ind <- sapply(uk_cases[,2], function(s) 
    !(s %in% c('Outside Wales','Unknown','...17','...18'))
)
uk_cases <- uk_cases[ind,]


# Assumes that metadata contains all areas in uk_cases and area names match perfectly
metadata <- read.csv("data/metadata.csv")
ind <- sapply(metadata[,1], function(s) 
    s %in% uk_cases[,2]
)
metadata <- metadata[ind,]


N <- nrow(uk_cases) # 149 number of regions in England & Wales only. 185 scotland too
D <- 100      # infection profile number of days
Tignore <- 3  # counts in most recent 3 days may not be reliable?
Tpred <- 1    # number of days held out for predictive probs eval
Tlik <- 7     # number of days for likelihood to infer Rt
Tall <- ncol(uk_cases)-2-Tignore  # number of days in time series used.
Tcond <- Tall-Tlik-Tpred       # number of days we condition on
Tproj <- 21              # number of days to project forward


Count <- uk_cases[1:N,3:(Tall+2)]

geoloc <- matrix(0, N, 2)
geodist <- matrix(0, N, N)
population <- rep(0.0, N)

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
    population[i] = metadata$POPULATION[j[1]]
  } else {
    print(sprintf("Cannot find region '%s'",region_name))
    for (r in 1:length(region_names)) {
      if (length(grep(region_names[r], region_name))>0) {
        geoloc[i, 1] = longitudes[r]
        geoloc[i, 2] = latitudes[r]
        population[i] = metadata$POPULATION[r]
        print(sprintf("...found region '%s'",region_names[r]))
      }
    }
  }
}

# compute distances between areas. Straight line, not actual travel distance. 
for (i in 1:N) {
  for (j in i:N) {
    # distance between two points on an ellipsoid (default is WGS84 ellipsoid), in units of 100km
    geodist[i, j] = distGeo(geoloc[i, 1:2], geoloc[j, 1:2]) / 100000
    geodist[j, i] = geodist[i, j]
  }
}

# compute fluxes for radiation model
flux <- matrix(0, N, N)
for (i in 1:N) {
  mi = population[i]
  js = order(geodist[i,])[2:N]
  nj = population[js]
  sij = cumsum(c(0.0,nj))[1:(N-1)]
  flux[i,js] = mi*nj/((mi+sij)*(mi+nj+sij))
  flux[i,i] = 1.0 - sum(flux[i,js]) # this statement weird, as fluxes don't sum to 1, unlike as claimed in Simini paper
}

Rmap_data <- list(
  N = N, 
  D = D, 
  Tall = Tall,
  Tcond = Tcond,
  Tlik = Tlik,
  Tproj = Tproj,
  Count = Count,
  geoloc = geoloc,
  geodist = geodist,
  flux = flux,
  infprofile = infprofile
  # local_sd = opt$local_sd,
  # global_sd = opt$global_sd,
  # gp_sd = opt$gp_sd,
  # gp_length_scale_sd = opt$gp_length_scale_sd
)

runname = sprintf('Rmap-%s-%s-%s-%s-%s', 
  opt$spatialkernel, 
  opt$localkernel, 
  opt$globalkernel, 
  opt$metapop, 
  opt$observation)
print(runname)


# copy the stan file and put in the right kernel
stan_file_name = paste('fits/', runname, '.stan', sep='')
content = readLines(paste('stan_files/', 'Rmap.stan',sep=''))
content = gsub(pattern="SPATIAL", replace=opt$spatialkernel, content)
content = gsub(pattern="LOCAL", replace=opt$localkernel, content)
content = gsub(pattern="GLOBAL", replace=opt$globalkernel, content)
content = gsub(pattern="METAPOP", replace=opt$metapop, content)
content = gsub(pattern="OBSERVATION", replace=opt$observation, content)
writeLines(content, stan_file_name)

fit <- stan(file = stan_file_name,
            data = Rmap_data, 
            iter = numiters, 
            chains = numchains,
            control = list(adapt_delta = .9))
# print(fit)

print(summary(fit, 
    pars=c("R0","gp_length_scale","gp_sigma","global_sigma","local_scale","precision","coupling_rate"), 
    probs=c(0.025, 0.25, 0.5, 0.75, 0.975))$summary)


s <- summary(fit, pars="Rt", probs=c(0.025, 0.25, .5, 0.75, .975))$summary
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
write.csv(df, paste('fits/', runname, '_logpred', '.csv', sep=''),
    row.names=FALSE)

inquotes <- function(s) paste('"',s,'"',sep='')

dates <- as.Date(colnames(uk_cases)[2+(Tcond+Tlik)], format='X%Y.%m.%d')
dates <- seq(dates,by=1,length.out=Tproj+1)[2:(Tproj+1)]
dates <- rep(dates,N)
areas <- rep(sapply(uk_cases[1:N,2],inquotes),Tproj)
dim(areas) <- c(N,Tproj)
areas <- t(areas)
dim(areas) <- c(N*Tproj)
df <- data.frame(area = areas, Date = dates, Cproj = Cproj)
colnames(df)[3:5] <- c("C_lower","C_median","C_upper")
df[,3:5] <- format(df[,3:5],digits=2)
write.csv(df, paste('fits/', runname, '_Cproj.csv', sep=''),
    row.names=FALSE,quote=FALSE)
write.csv(df, paste('website/Cproj.csv', sep=''),
    row.names=FALSE,quote=FALSE)

areas <- sapply(uk_cases[1:N,2],inquotes)
df <- data.frame(area = areas, Rt = Rt)
colnames(df)[2:4] <- c("Rt_lower","Rt_median","Rt_upper")
df[,2:4] <- format(df[,2:4],digits=2)
write.csv(df, paste('fits/', runname, '_Rt.csv', sep=''),
    row.names=FALSE,quote=FALSE)
write.csv(df, paste('website/Rt.csv', sep=''),
    row.names=FALSE,quote=FALSE)



# df <- data.frame(area = uk_cases[1:N,2], Rt = Rt, Cproj = Cproj)
# colnames(df) <- c("area","Rtlower","Rtmedian","Rtupper","Cprojlower","Cprojmedian","Cprojupper")
# write.csv(df, paste('fits/', runname, '_RtCproj', '.csv', sep=''),row.names=FALSE)

saveRDS(fit, paste('fits/', runname, '_stanfit', '.rds', sep=''))

print(runname)

pairs(fit, pars=c("R0","gp_length_scale","gp_sigma","global_sigma","local_scale","precision","coupling_rate","rad_prob"))

