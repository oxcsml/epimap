Rmap_merge = function(env,cleaned_sample_ids) { 
env$cleaned_sample_ids = cleaned_sample_ids
with(env, {

writemergedresults = function(data,filename,...) {
  write.csv(data,sprintf('%s/merged_%s.csv',opt$results_directory,filename),...)
}

numruns = length(cleaned_sample_ids)
load_samples = function(pars) {
  samples = do.call(rbind, lapply(1:numruns, function(i) {
    read.csv(paste(
      opt$results_directory,
      '/',cleaned_sample_ids[i],
      '_',pars,
      '_samples.csv',
      sep=''
    ))
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

Rt = Rt[sapply(1:N,function(i)rep((i-1)*Mstep+c(1:Mstep,rep(Mstep,Mproj)),each=Tstep)),]
df <- area_date_dataframe(
    quoted_areas, 
    days_all,
    provenance,
    format(round(Rt,2),nsmall=2),
    #c("Rt_10","Rt_20","Rt_30","Rt_40","Rt_50","Rt_60","Rt_70","Rt_80","Rt_90")
    c("Rt_2_5","Rt_10","Rt_20","Rt_25","Rt_30","Rt_40","Rt_50",
      "Rt_60","Rt_70","Rt_75","Rt_80","Rt_90","Rt_97_5")
)
writemergedresults(df, 'Rt', row.names=FALSE, quote=FALSE)

#################################################################
# Rt exceedance probabilities
thresholds = c(.8, .9, 1.0, 1.1, 1.2, 1.5, 2.0)
numthresholds = length(thresholds)
numsamples = numruns * numiters/2
Rt <- as.matrix(Rt_samples)
dim(Rt) <- c(numsamples,Mstep,N)
Pexceedance = array(0.0,dim=c(Mstep,N,numthresholds))
for (k in 1:Mstep) {
  for (i in 1:N) {
    for (x in 1:numthresholds) {
      Pexceedance[k,i,x] = mean(Rt[,k,i]>thresholds[x])
    }
  }
}
Pexceedance = Pexceedance[c(1:Mstep,rep(Mstep,Mproj)),,]
Pexceedance <- Pexceedance[sapply(1:(Mstep+Mproj),function(k)rep(k,Tstep)),,]
dim(Pexceedance) <- c(Tstep*(Mstep+Mproj)*N,numthresholds)
df <- area_date_dataframe(
    quoted_areas, 
    days_all,
    provenance,
    format(round(Pexceedance,2),nsmall=2),
    c("P_08","P_09","P_10","P_11","P_12","P_15","P_20")
)
writemergedresults(df, 'Pexceed', row.names=FALSE, quote=FALSE)



#################################################################
# posterior predictives and projections
Cpred = t(apply(Cpred_samples,2,quantile,
    probs=c(0.025, 0.25, .5, 0.75, .975)
))

df <- area_date_dataframe(
    quoted_areas,
    seq(dates[Tcond]+1,by=1,length.out=Mstep*Tstep),
    rep('inferred',Tlik),
    format(round(Cpred,1),nsmall=1),
    #c("C_2_5","C_25","C_50","C_75","C_97_5")
    c("C_025","C_25","C_50","C_75","C_975")
)
writemergedresults(df, 'Cpred', row.names=FALSE, quote=FALSE)

####################################################################################
# weekly counts. Includes 1 last column of actual counts among days ignored in model
Cweekly <- as.matrix(Count[,(Tcond+1):(Tcond+Tlik)])
dim(Cweekly) <- c(N,Tstep,Mstep)
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
Cweekly <- Cweekly[sapply(1:(Mstep+Mproj),function(k)rep(k,Tstep)),]
dim(Cweekly) <- c(N*(Tlik+Tstep*Mproj))
df <- area_date_dataframe(
    quoted_areas,
    days_all,
    provenance,
    format(Cweekly,digits=3),
    c("C_weekly")
)
writemergedresults(df, 'Cweekly', row.names=FALSE, quote=FALSE)



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
writemergedresults(df, 'Cproj', row.names=FALSE, quote=FALSE)

})
}
