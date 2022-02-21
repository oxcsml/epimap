
# compute fluxes for radiation model.
# relax radiation model by adding noise to distance with sd 10km, and bootstrapping

compute_radiation_fluxes = function(numrep, radiation_length_scale, population, geodist) {
  flux <- array(0, dim=c(N, N))
  flux_var <-  array(0, dim=c(N, N))
  for (r in 1:numrep) {
    ff <- matrix(0, N, N)
    for (i in 1:N) {
      mi = population[i]
      if (i==1) {
        jj <- 2:N
      } else if (i==N) {
        jj <- 1:(N-1)
      } else {
        jj = c(1:(i-1),(i+1):N)
      }
      js = jj[order(geodist[i,jj]+rnorm(N-1,sd=radiation_length_scale))]
      nj = population[js]
      sij = cumsum(c(0.0,nj))[1:(N-1)]
      ff[i,js] = mi*nj/((mi+sij)*(mi+nj+sij))
      ff[i,i] = 1.0 - sum(ff[i,js]) # this statement weird, as fluxes don't sum to 1, unlike as claimed in Simini paper
    }
    flux = flux + ff
    flux_var = flux_var + ff^2
  }
  flux = flux / numrep
  flux_sd = (flux_var - numrep*(flux^2))/(numrep-1)
  flux_sd[flux_sd<0]=0
  flux_sd = sqrt(flux_sd)/sqrt(numrep)

  list("flux"=flux, "flux_sd"=flux_sd)
}