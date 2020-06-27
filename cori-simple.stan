// saved as schools.stan
data {
  int<lower=1> N;           // number of regions
  int<lower=1> T;           // number of days
  int<lower=1> T0;          // number of recent days for likelihood computation
  int<lower=1> D;           // length of infection profile. 
  int C[N, T];              // case counts
  vector[D] infprofile;     // infection profile
}

transformed data {
  vector[D] infprofile_rev; // reversed infection profile
  vector[T0] convout[N];       // convolution between C and infprofile
  vector[T] Cvec[N];
  for (i in 1:D)
    infprofile_rev[i] = infprofile[D-i+1];
  for (j in 1:N) {
    for (i in 1:T)
      Cvec[j][i] = C[j,i];
    for (i in 1:T0) {
      int L = min(D,T-T0+i-1); // length of infection profile that overlaps with case counts in past days
      convout[j][i] = dot_product(Cvec[j][T-T0-L+i:T-T0+i-1], infprofile_rev[D-L+1:D]);
    }
  }
}

parameters {
  real<lower=0> Rt[N];                // instantaneous reproduction number
}

model {
  for (j in 1:N) {
    for (i in 1:T0) {
      C[j][T-T0+i] ~ poisson(Rt[j] * convout[j][i]);
    }
  }
}

