functions {
  matrix exp_quad(matrix dist, real func_sigma, real length_scale) {
    return square(func_sigma) * exp( (-1.0 / (2.0 * square(length_scale))) * square(dist));
  }
  matrix matern12(matrix dist, real func_sigma, real length_scale) {
    return square(func_sigma) * exp(- dist / length_scale);
  }
  matrix matern32(matrix dist, real func_sigma, real length_scale) {
    return square(func_sigma) * (1 + ((sqrt(3) * dist) / length_scale)) .* exp(- (sqrt(3) * dist) / length_scale);
  }
  matrix matern52(matrix dist, real func_sigma, real length_scale) {
    return square(func_sigma) * (1 + ((sqrt(5) * dist) / length_scale) + ((5 * dist .* dist)/(3 * square(length_scale)))) .* exp(- (sqrt(5) * dist) / length_scale);
  }
}

data {
  int<lower=1> N;           // number of regions
  int<lower=1> D;           // length of infection profile
  int<lower=1> Tall;        // number of all days in case count time series
  int<lower=1> Tcond;       // number of days we will condition on
  int<lower=1> Tlik;        // number of days for likelihood computation
  int<lower=0> Tproj;       // number of days to forecast
  int Count[N, Tall];       // case counts
  vector[2] geoloc[N];      // geo locations of regions
  matrix[N,N] geodist;      // geo locations of regions
  vector[D] infprofile;     // infection profile aka serial interval distribution
}

transformed data {
  int Tcur = Tcond+Tlik;    // index of day on which we are estimating Rt
  int Tpred = Tall-Tcur;    // number of days to calculate predictive probabilities for

  vector[Tall] Creal[N];     // real type version of Count
  vector[D] infprofile_rev; // reversed infection profile

  // precompute convolutions between Count and infprofile 
  matrix[N,Tlik] convlik;      // for use in likelihood computation
  matrix[N,Tpred] convpred;    // for use in predictive probs of future counts
  matrix[N,Tproj] convproj;    // for use in forecasting into future 

  // reverse infection profile
  for (i in 1:D)
    infprofile_rev[i] = infprofile[D-i+1];

  for (j in 1:N) {
    for (i in 1:Tall)
      Creal[j,i] = Count[j,i];

    // precompute convolutions between counts and infprofile
    for (i in 1:Tlik) {
      int L = min(D,Tcond+i-1); // length of infection profile that overlaps with case counts 
      convlik[j,i] = dot_product(Creal[j][Tcond-L+i:Tcond-1+i], infprofile_rev[D-L+1:D]);
    }
    for (i in 1:Tpred) {
      int L = min(D,Tcur+i-1); // length of infection profile that overlaps with case counts 
      convpred[j,i] = dot_product(Creal[j][Tcur-L+i:Tcur-1+i], infprofile_rev[D-L+1:D]);
    }
    for (i in 1:Tproj) {
      int L = min(D,Tcur+i-1); // length of infection profile that overlaps with case counts 
      convproj[j,i] = dot_product(Creal[j][Tcur-L+i:Tcur], infprofile_rev[D-L+1:D-i+1]);
    }
  }
}

parameters {
  real<lower=0> length_scale;
  real<lower=0> func_sigma;
  real<lower=0> data_sigma;
  // real<lower=0> avg_sigma;
  vector[N] eta;

  real<lower=0> dispersion;
  real<lower=0> Ravg;
  real<lower=0,upper=1> immigration_rate;
}

transformed parameters {
  vector[N] Rt;                 // instantaneous reproduction number

  {
    real data_sigma2 = square(data_sigma);
    // real avg_sigma2 = square(avg_sigma);
    matrix[N,N] K;
    matrix[N,N] L;

    // GP prior.

    // K = cov_exp_quad(geoloc, func_sigma, length_scale); // kernel
    // K = KERNEL(geodist, func_sigma, length_scale); // kernel
    K = square(func_sigma) * exp( (-1.0 / (2.0 * square(length_scale))) * square(geodist));
    for (i in 1:N) {
      K[i,i] = K[i,i] + data_sigma2;
      // for (j in 1:N) {
      //   K[i,j] = K[i,j] + avg_sigma2;
      // }
    }

    L = cholesky_decompose(K);
    // Rt = exp(L * eta);
    Rt = Ravg * exp(L * eta);
  }
}

model {
  vector[Tlik] immigration;

  Ravg ~ normal(1.0,1.0);
  immigration_rate ~ normal(.05, .05);
  dispersion ~ normal(0,5);

  // GP prior density
  length_scale ~ normal(0.0,10.0);
  func_sigma ~ normal(0.0, 1.0);
  data_sigma ~ normal(0.0, 1.0);
  // avg_sigma ~ normal(0.0, 1.0);
  eta ~ std_normal();


  for (i in 1:Tlik) {
    immigration[i] = 0.0;
    for (j in 1:N)
      immigration[i] += Rt[j] * convlik[j][i];
    immigration[i] = immigration[i] / N;
  }

  for (j in 1:N) {
    for (i in 1:Tlik) {
      Count[j,Tcond+i] ~ neg_binomial_2(
          (1.0-immigration_rate) *  Rt[j] * convlik[j,i] +
          immigration_rate * immigration[i], dispersion);
    }
  }

/*
  for (i in 1:Tlik) {
    real immigration; // cross-regional infections (immigration)
    immigration = 0.0;
    for (j in 1:N) 
      immigration += Rt[j] * convlik[j,i]; 
    immigration = immigration / N;
   
    // likelihood for inferring Rt 
    for (j in 1:N) {
      Count[j,Tcond+i] ~ neg_binomial_2(
          (1.0-immigration_rate) *  Rt[j] * convlik[j,i] +
          immigration_rate * immigration, 
          dispersion);
    }
  }
*/
}

generated quantities {
  // real Ravg = mean(Rt);
  vector[Tpred] Ppred[N];
  vector[Tproj] Cproj[N]; // faster than matrix

  // predictive probability of future counts
  for (i in 1:Tpred) {
    real immigration;
    immigration = 0.0;
    for (j in 1:N) 
      immigration += Rt[j] * convpred[j,i]; 
    immigration = immigration / N;

    for (j in 1:N) {
      Ppred[j,i] = exp(neg_binomial_2_lpmf(Count[j,Tcur+i] |
          (1.0-immigration_rate) *  Rt[j] * convpred[j,i] +
          immigration_rate * immigration, 
          dispersion));
    }
  }
  // forecasting *mean* counts given parameters
  for (i in 1:Tproj) {
    real immigration;
    vector[N] convprojall;
    immigration = 0.0;
    for (j in 1:N) {
      convprojall[j] = convproj[j,i] + 
          dot_product(Cproj[j][1:(i-1)], infprofile_rev[D-i+2:D]);
      immigration += Rt[j] * convprojall[j];
    }
    immigration = immigration / N;

    for (j in 1:N) {
      // Cproj[j,i] = neg_binomial_2_rng( 
      //     (1.0-immigration_rate) *  Rt[j] * convprojall[j] +
      //     immigration_rate * immigration, 
      //     dispersion);
      Cproj[j,i] = ( 
          (1.0-immigration_rate) *  Rt[j] * convprojall[j] +
          immigration_rate * immigration);
    }

  }

}


