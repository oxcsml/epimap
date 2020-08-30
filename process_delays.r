
library(rstan)

numchains = 1
numiters = 3000

options(mc.cores = min(numchains,parallel::detectCores()))
rstan_options(auto_write = TRUE)

source('read_data.r')
Tignore <- 5  # counts in most recent 7 days may not be reliable?
Tall <- Tall-Tignore
Count <- Count[,1:Tall]

shortlist <- c(
  'Herefordshire, County of',
  'Oldham',
  'Birmingham',
  'Manchester',
  'Northampton',
  'Slough',
  'Torridge', #7
  'Teignbridge',
  'Highland'
)
area <- shortlist[7]
print(area)

Nstep <- 15
Tstep <- 7

Tip <- 30
infprofile <- infprofile[1:Tip]
infprofile <- infprofile/sum(infprofile)

Tdp <- 14
Adp <- 5.0
Bdp <- 1.0
delayprofile <- pgamma(1:Tdp,shape=Adp,rate=Bdp)
delayprofile <- delayprofile/delayprofile[Tdp]
delayprofile <- delayprofile - c(0.0,delayprofile[1:(Tdp-1)])

#############################################################################################
# Main computation

Rmap_local_smooth_data <- list(
  Tall = Tall,
  Tstep = Tstep, 
  Nstep = Nstep,
  Count = Count[area,],
  Tip = Tip,
  infprofile = infprofile,
  Tdp = Tdp,
  delayprofile = delayprofile,
  mu_scale = 0.5,
  sigma_scale = 0.5,
  alpha_scale = 0.25,
  phi_latent_scale = 10.0,
  phi_observed_scale = 10.0
)

start_time <- Sys.time()

fit <- stan(file = 'stan_files/Rmap-local-smooth.stan',
            data = Rmap_local_smooth_data, 
            iter = numiters, 
            chains = numchains,
            control = list(adapt_delta = .9))

end_time <- Sys.time()

print("Time to run")
print(end_time - start_time)


#################################################################
# Summary of fit
print(
    summary(fit, 
        pars=c("mu","sigma","alpha","phi_latent","phi_observed","Rt"),
        probs=c(.025,0.5,.975)
    )$summary
)
 
####################################################################
# pairs plot
pdf(paste("local-smoothed/pairs-",area,".pdf",sep=""),width=9,height=9)
pairs(fit, pars=c("mu","sigma","alpha","phi_latent","phi_observed"))
dev.off()

pdf(paste("local-smoothed/samples-",area,".pdf",sep=""),width=9,height=9)
par(mfrow=c(4,2))
par(oma=c(0,0,0,0))
par(mar=c(1,1,1,1))
CinferCI = summary(fit,pars="Cinfer",probs=c(0.025,0.25,0.5,0.75,0.975))$summary
ind = (Tall-Nstep*Tstep+1):Tall
CinferCI = CinferCI[,c("2.5%","50%","97.5%")]
Cinfer = extract(fit,pars="Cinfer",permuted=FALSE)
dim(Cinfer) <- c(numiters/2,Nstep*Tstep)
for (i in seq(from=187,by=187,to=numiters/2)) {
  #dev.new()
  plot(t(Count[area,]),pch=20,ylim=c(0,max(Count[area,ind])))
  for (j in 1:3) {
    lines(ind,CinferCI[,j])
  }
  points(ind,Cinfer[i,],col='red',pch=20)
  #dev.off()
  #a <- readline()
}
dev.off()

