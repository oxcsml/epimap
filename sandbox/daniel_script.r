source('epimap/epiclean.r')

m = c(rep(.5,21*2)^2,1.3^seq(1,to=21),1.3^21*0.9^seq(1,to=21),rep(1.3^21*.9^21,21))
Count = rpois(21*5,m)
print(Count)

Tcond = 21
Nstep = 21*4
Tstep = 1
Nproj = 21

infprofile = read.csv('data/serial_interval.csv')$fit
Tip = 30
infprofile <- infprofile[1:Tip]
infprofile <- infprofile/sum(infprofile)

Tdp <- 21
presymptomdays <- 2
Tdpnz <- Tdp - presymptomdays
Adp <- 5.8
Bdp <- 0.948
testdelayprofile <- pgamma(1:Tdpnz,shape=Adp,rate=Bdp)
testdelayprofile <- testdelayprofile/testdelayprofile[Tdpnz]
testdelayprofile <- testdelayprofile - c(0.0,testdelayprofile[1:(Tdpnz-1)])
testdelayprofile <- c(rep(0, presymptomdays), testdelayprofile)

gp_time_scale = 14
gp_time_decay_scale = .5
num_iterations = 3000
fit = epiclean(
  Count = Count,
  Tcond = Tcond,
  Nstep = Nstep,
  Nproj = Nproj,
  Tstep = Tstep,
  infprofile = infprofile,
  testdelayprofile = testdelayprofile,
  gp_time_scale = gp_time_scale,
  gp_time_decay_scale = gp_time_decay_scale,
  num_iterations = num_iterations
)
print(summary(fit$stanfit,
  pars=c(
    "mu",
    "sigma",
    "alpha",
    "gp_time_length_scale",
    "phi_latent",
    "phi_observed",
    "xi",
    "Noutliers",
    "meandelay"),
  probs=c(.025,.5,.975)
)$summary)

print(fit$Rt[,c("2.5%","50%","97.5%")])

