// saved as schools.stan
data {
  int<lower=1> N;           // number of regions
  int<lower=1> T;           // number of days
  int<lower=1> T0;          // number of recent days for likelihood computation
  int<lower=1> D;           // length of infection profile
  int<lower=0> Tproj;       // number of days to forward project
  int C[N, T];              // case counts
  vector[2] geoloc[N];       // geo locations of regions
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
  real<lower=0> length_scale;
  real<lower=0> func_sigma;
  real<lower=0> data_sigma;
  vector[N] eta;

  real<lower=0> dispersion;
  real<lower=0> Ravg;
  real<lower=0,upper=1> immigration_rate;
}

transformed parameters {
  vector[N] Rt;                 // instantaneous reproduction number

  {
    real data_sigma2;
    matrix[N,N] K;
    matrix[N,N] L;

    data_sigma2 = square(data_sigma);

    K = cov_exp_quad(geoloc, func_sigma, length_scale); // kernel
    for (i in 1:N)
      K[i,i] = K[i,i] + data_sigma2;

    L = cholesky_decompose(K);
    Rt = Ravg * exp(L * eta);
  }
}

model {
  vector[T0] immigration;

  Ravg ~ normal(1.0,1.0);
  immigration_rate ~ beta(1.0 * .01, 1.0 * .99);
  dispersion ~ normal(0,5);

  // GP
  length_scale ~ normal(0.1,1.0);
  func_sigma ~ normal(0.1, 1.0);
  data_sigma ~ normal(0.1, 1.0);
  eta ~ std_normal();

  for (i in 1:T0) {
    immigration[i] = 0.0;
    for (j in 1:N) 
      immigration[i] += Rt[j] * convout[j][i]; 
    immigration[i] = immigration[i] / N;
  }
  
  for (j in 1:N) {
    for (i in 1:T0) {
      C[j][T-T0+i] ~ neg_binomial_2(
          (1.0-immigration_rate) *  Rt[j] * convout[j][i] +
          immigration_rate * immigration[i], dispersion);
    }
  }
}

generated quantities {
  vector[Tproj] Cproj[N];

  for (j in 1:N) {
    for (i in 1:Tproj) {
      int L = min(D,T+i-1); // length of infection profile that overlaps with case counts in past days
      Cproj[j][i] = Rt[j] * (dot_product(Cvec[j][T+i-L:T], infprofile_rev[D-L+1:D-i+1]) + 
                            dot_product(Cproj[j][1:i-1], infprofile_rev[D-i+2:D]));
    }
  }

}

