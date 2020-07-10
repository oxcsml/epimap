library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

infprofile <- read.csv("serial_interval.csv")$fit

uk_cases <- read.csv("uk_cases.csv")

metadata <- read.csv("metadata.csv")

T <- ncol(uk_cases)-2-7  # number of days; last 7 days counts not reliable
N <- 149                 # number of regions in England only
T0 <- 7                  # number of days to average over to estimate Rt
D <- 100                 # infection profile number of days
Tproj <- 21              # number of days to project forward


C <- uk_cases[1:N,3:(T+2)]

geoloc <- matrix(0, N, 2)

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
  

cori_dat <- list(N = N, T = T, T0 = T0, Tproj = Tproj, D = D, C = C,  
                 geoloc = geoloc,
                 infprofile = infprofile)

# fit <- stan(file = 'cori-simple.stan', data = cori_dat)
# fit <- stan(file = 'cori-gp.stan', data = cori_dat)
fit <- stan(file = 'cori-gp-immi.stan', data = cori_dat, iter=4000, control = list(adapt_delta = .9))
# print(fit)

print(summary(fit, pars=c("Ravg","length_scale","func_sigma","data_sigma","dispersion","immigration_rate"), probs=c(0.025, 0.5, 0.975))$summary)


s <- summary(fit, pars="Rt", probs=c(0.025, .5, .975))$summary
Rt <- s[,c("2.5%","50%","97.5%")]
Rt <- t(t(Rt))

sprintf("median Rt range: [%f, %f]",min(Rt[,2]),max(Rt[,2]))

s <- summary(fit, pars="Cproj", probs=c(0.025, .5, .975))$summary
Cproj <- s[,c("2.5%","50%","97.5%")]
Cproj <- t(t(Cproj))
Cproj <- Cproj[Tproj*(1:N),]

sprintf("median Cproj range: [%f, %f]",min(Cproj[,2]),max(Cproj[,2]))

df <- data.frame(area = uk_cases[1:N,2], Rt = Rt, Cproj = Cproj)
colnames(df) <- c("area","Rtlower","Rtmedian","Rtupper","Cprojlower","Cprojmedian","Cprojupper")

write.csv(df,"RtCproj.csv",row.names=FALSE)

