
library(rstan)

numchains = 1
numiters = 3000

options(mc.cores = min(numchains,parallel::detectCores()))
rstan_options(auto_write = TRUE)

source('read_data.r')
Tignore <- 5  # counts in most recent 5 days may not be reliable?
Tall <- Tall-Tignore
Count <- Count[,1:Tall]

Nsample <- 8

Tip <- 30
infprofile <- infprofile[1:Tip]
infprofile <- infprofile/sum(infprofile)

Tdp <- 14
Adp <- 5.0
Bdp <- 1.0
delayprofile <- pgamma(1:Tdp,shape=Adp,rate=Bdp)
delayprofile <- delayprofile/delayprofile[Tdp]
delayprofile <- delayprofile - c(0.0,delayprofile[1:(Tdp-1)])

Tstep <- 7
Nstep <- floor((Tall-max(Tip,Tdp)) / Tstep)
Tlik <- Nstep*Tstep
Tcond <- Tall-Tlik


#############################################################################################

shortlist <- c(
  'Orkney',
  'Herefordshire, County of',
  'Oldham',
  'Birmingham',
  'Manchester',
  'Northampton',
  'Slough',
  'Torridge', 
  'Teignbridge',
  'Highland'
)
do_shortlist = FALSE
if (do_shortlist) {
  quoted_areas <- shortlist
  areas <- shortlist
  Count <- Count[shortlist,]
  N = length(shortlist)
}

Cinfer_sample <- array(0,c(N,Tall,Nsample))
Cinfer_mean <- array(0,c(N,Tall))
Clean_sample <- array(0,c(N,Tall,Nsample))
Cinfer_mean[,1:Tcond] = as.matrix(Count[,1:Tcond])
for (i in 1:Nsample) {
  Cinfer_sample[,1:Tcond,i] = as.matrix(Count[,1:Tcond])
  Clean_sample[,1:Tcond,i] = as.matrix(Count[,1:Tcond])
}

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
  phi_observed_scale = 10.0,
  outlier_threshold = .95,
  exogeneous_infections = 1e-4
)

init = list()
init[[1]] = list(
  mu = 0.0,
  sigma = 0.01,
  alpha1 = .95,
  phi_latent = 1.0,
  phi_observed = 5.0,
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
        pars=c("mu","sigma","alpha","phi_latent","phi_observed","Noutliers","Rt"),
        probs=c(0.5)
    )$summary
)
 
skip = numiters/2/Nsample;
####################################################################
Cinfer_s <- extract(fit,pars="Cinfer",permuted=FALSE)
Cinfer_s <- Cinfer_s[seq(skip,by=skip,length.out=Nsample),,];
dim(Cinfer_s) <- c(Nsample,Tlik)
Cinfer_s <- t(Cinfer_s)
dim(Cinfer_s) <- c(1,Tlik,Nsample)
Cinfer_sample[area_index,(Tcond+1):(Tcond+Tlik),] <- Cinfer_s
Cinfer_m <- summary(fit,pars="Cinfer",probs=c(0.5))$summary
Cinfer_m <- t(as.matrix(Cinfer_m[,"mean"]))
Cinfer_mean[area_index,(Tcond+1):(Tcond+Tlik)] <- Cinfer_m

####################################################################
Clean_s <- extract(fit,pars="Clean",permuted=FALSE)
Clean_s <- Clean_s[seq(skip,by=skip,length.out=Nsample),,]
Clean_s <- t(Clean_s)
dim(Clean_s) <- c(1,Tall,Nsample)
Clean_sample[area_index,,] <- Clean_s

####################################################################
# pairs plot
pdf(paste("local-cleaned/pairs-",area,".pdf",sep=""),width=9,height=9)
pairs(fit, pars=c("mu","sigma","alpha","phi_latent","phi_observed"))
dev.off()

pdf(paste("local-cleaned/Cinfer-",area,".pdf",sep=""),width=9,height=9)
par(mfrow=c(4,2))
par(oma=c(0,0,0,0))
par(mar=c(1,1,1,1))
CinferCI = summary(fit,pars="Cinfer",probs=c(0.025,0.25,0.5,0.75,0.975))$summary
ind = (Tall-Nstep*Tstep+1):Tall
CinferCI = CinferCI[,c("2.5%","50%","97.5%")]
for (i in 1:Nsample) {
  plot(t(Count[area,]),pch=20,ylim=c(0,max(Count[area,ind])))
  for (j in 1:3) {
    lines(ind,CinferCI[,j])
  }
  points(ind,Cinfer_sample[area_index,ind,i],col='red',pch=20)
}
dev.off()

pdf(paste("local-cleaned/clean-",area,".pdf",sep=""),width=9,height=9)
par(mfrow=c(4,2))
par(oma=c(0,0,0,0))
par(mar=c(1,1,1,1))
CleanCI = summary(fit,pars="Clean",probs=c(0.025,0.25,0.5,0.75,0.975))$summary
ind = (Tall-Nstep*Tstep+1):Tall
CleanCI = CleanCI[,c("2.5%","50%","97.5%")]
for (i in 1:Nsample) {
  plot(t(Count[area,]),pch=20,ylim=c(0,max(Count[area,ind])))
  for (j in 1:3) {
    lines(ind,CleanCI[ind,j])
  }
  points(ind,Clean_sample[area_index,ind,i],col='red',pch=20)
}
dev.off()


}

rownames(Cinfer_mean) <- quoted_areas
write.csv(Cinfer_mean,'data/Cinfer_mean.csv',quote=FALSE)
for (i in 1:Nsample) {
  cc <- Cinfer_sample[,,i]
  rownames(cc) <- quoted_areas
  colnames(cc) <- colnames(Count)
  write.csv(cc,paste('data/Cinfer_sample',i,'.csv',sep=''),quote=FALSE)
  cc <- Clean_sample[,,i]
  rownames(cc) <- quoted_areas
  colnames(cc) <- colnames(Count)
  write.csv(cc,paste('data/Clean_sample',i,'.csv',sep=''),quote=FALSE)
}
