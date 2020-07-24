
kernel = c('none','none','matern12','matern32','exp_quad','none','none','matern12','matern12')
metapop = c('none','none','none','none','none','uniform1','uniform2','uniform1','uniform1')
obs = c('poisson','negative_binomial','negative_binomial','negative_binomial','negative_binomial','negative_binomial','negative_binomial','negative_binomial','negative_binomial')

minval = -log(200)

dbaseline = pmax(minval,data.matrix(read.csv(sprintf("fits/cori_local_logpred.csv"))[,2:8]))
n = length(dbaseline)

for (i in 1:9) {
  d = pmax(minval,data.matrix(read.csv(sprintf("fits/Rmap-%s-%s-%s_logpred.csv",kernel[i],metapop[i],obs[i]))[,2:8]))
  diff = d-dbaseline

  m = mean(diff)
  me = qt(.95,n-1) * sd(diff) / sqrt(n)
  print(sprintf("%8s %8s %17s: %f, CI[%f, %f]",kernel[i],metapop[i],obs[i],m,m-me,m+me))
}

