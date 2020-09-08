
library(rstan)

numchains = 1
numiters = 3000

options(mc.cores = min(numchains,parallel::detectCores()))
rstan_options(auto_write = TRUE)

source('read_data.r')
Tignore <- 5  # counts in most recent 7 days may not be reliable?
Tall <- Tall-Tignore
Count <- Count[,1:Tall]

Nstep <- 15
Tstep <- 7
Tlik <- Nstep*Tstep
Tcond <- Tall-Tlik

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

Cinfer_sample <- Count
Cinfer_mean <- Count

for (area_index in 1:N) {

area = areas[area_index]
print(area)


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

init = list()
init[[1]] = list(
  mu = 0,
  sigma = 0.1,
  alpha1 = .95,
  phi_latent = 2,
  phi_observed = 2,
  'Reta[1]' = 1,
  'Reta[2]' = 1,
  'Reta[3]' = 1,
  'Reta[4]' = 1,
  'Reta[5]' = 1,
  'Reta[6]' = 1,
  'Reta[7]' = 1,
  'Reta[8]' = 1,
  'Reta[9]' = 1,
  'Reta[10]' = 1,
  'Reta[11]' = 1,
  'Reta[12]' = 1,
  'Reta[13]' = 1,
  'Reta[14]' = 1,
  'Reta[15]' = 1
)

start_time <- Sys.time()

fit <- stan(file = 'stan_files/Rmap-local-smooth.stan',
            data = Rmap_local_smooth_data, 
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
        pars=c("mu","sigma","alpha","phi_latent","phi_observed","Rt"),
        probs=c(0.5)
    )$summary
)
 
####################################################################
Cinfer_s <- as.matrix(extract(fit,pars="Cinfer",permuted=TRUE))[[1]]
Cinfer_sample[area,(Tcond+1):(Tcond+Tlik)] <- Cinfer_s[1,]
Cinfer_m <- summary(fit,pars="Cinfer",probs=c(0.5))$summary
Cinfer_m <- Cinfer_m[,"mean"]
Cinfer_mean[area,(Tcond+1):(Tcond+Tlik)] <- Cinfer_m

####################################################################
# pairs plot
pdf(paste("local-cleaned/pairs-",area,".pdf",sep=""),width=9,height=9)
pairs(fit, pars=c("mu","sigma","alpha","phi_latent","phi_observed"))
dev.off()

pdf(paste("local-cleaned/samples-",area,".pdf",sep=""),width=9,height=9)
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

}

rownames(Cinfer_sample) <- quoted_areas
write.csv(Cinfer_sample,'data/Clean_sample.csv',sep=' ',quote=FALSE)
rownames(Cinfer_mean) <- quoted_areas
write.csv(Cinfer_mean,'data/Clean_mean.csv',sep=' ',quote=FALSE)

