functions {
  // GP kernel choices

  matrix exp_quad_kernel(matrix dist, real gp_sigma, real gp_length_scale) {
    return square(gp_sigma) * exp( (-1.0 / (2.0 * square(gp_length_scale))) * square(dist));
  }
  matrix matern12_kernel(matrix dist, real gp_sigma, real gp_length_scale) {
    return square(gp_sigma) * exp(- dist / gp_length_scale);
  }
  matrix matern32_kernel(matrix dist, real gp_sigma, real gp_length_scale) {
    return (square(gp_sigma) + 
        (square(gp_sigma)*sqrt(3.0)/gp_length_scale) * dist) .* 
        exp((-sqrt(3.0)/gp_length_scale) * dist) ;
  }
  matrix matern52_kernel(matrix dist, real gp_sigma, real gp_length_scale) {
    return square(gp_sigma) * 
        (1.0 + ((sqrt(5.0)/gp_length_scale) * dist)  + 
        ((5.0/(3.0 * square(gp_length_scale))) * dist .* dist)) .* 
        exp(- (sqrt(5.0)/gp_length_scale)* dist);
  }

  matrix none_kernel(matrix dist, real gp_sigma, real gp_length_scale) {
    return rep_matrix(0.0, rows(dist), cols(dist));
  }

  // Meta-population infection rate model choices

  matrix uniform1_metapop(
      vector Rt, matrix convlik, real coupling_rate) {
    int T = cols(convlik);
    int N = rows(convlik);
    row_vector[T] convavg;
    matrix[N,T] convout;
    
    for (i in 1:T) {
      convavg[i] = 0.0;
      for (j in 1:N)
        convavg[i] += Rt[j] * convlik[j,i];
      convavg[i] = convavg[i] / N;
    }

    for (j in 1:N)
      for (i in 1:T)
        convout[j,i] = (
            (1.0-coupling_rate) *  Rt[j] * convlik[j,i] +
            coupling_rate * convavg[i]
        );
    return convout;
  }

  matrix uniform2_metapop(
      vector Rt, matrix convlik, real coupling_rate) {
    int T = cols(convlik);
    int N = rows(convlik);
    row_vector[T] convavg;
    matrix[N,T] convout;
    
    for (i in 1:T) {
      convavg[i] = 0.0;
      for (j in 1:N)
        convavg[i] += convlik[j,i];
      convavg[i] = convavg[i] / N;
    }

    for (j in 1:N)
      for (i in 1:T)
        convout[j,i] = Rt[j] * (
            (1.0-coupling_rate) * convlik[j,i] +
            coupling_rate * convavg[i]
        );
    return convout;
  }

  matrix none_metapop(
      vector Rt, matrix convlik, real coupling_rate) {
    int T = cols(convlik);
    int N = rows(convlik);
    matrix[N,T] convout;

    for (j in 1:N)
      for (i in 1:T)
        convout[j,i] = Rt[j] * convlik[j,i];

    return convout;
  }

  // Case count likelihood choices

  real poisson_likelihood_lpmf(int count, real mean, real dispersion) {
    return poisson_lpmf(count | mean);
  }

  real negative_binomial_likelihood_lpmf(int count, real mean, real dispersion) {
    return neg_binomial_2_lpmf(count | mean, dispersion);
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
  vector[D] infprofile;     // infection profile aka serial interval distribution
  matrix[N,N] geodist;      // distance between locations
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
  real<lower=0> gp_length_scale;
  real<lower=0> gp_sigma;
  real<lower=0> local_sigma;
  real<lower=0> global_sigma;
  vector[N] eta;

  real<lower=0> dispersion;
  // real<lower=0> Ravg;
  real<lower=0,upper=1> coupling_rate;
}

transformed parameters {
  vector[N] Rt;                 // instantaneous reproduction number

  {
    real local_sigma2 = square(local_sigma);
    real global_sigma2 = square(global_sigma);
    matrix[N,N] K;
    matrix[N,N] L;

    // GP prior.

    K = KERNEL_kernel(geodist, gp_sigma, gp_length_scale); // kernel
    for (i in 1:N) {
      K[i,i] = K[i,i] + local_sigma2;
    }
    K = K + global_sigma2;

    L = cholesky_decompose(K);
    Rt = exp(L * eta);
    // Rt = Ravg * exp(L * eta);
  }
}

model {
  vector[Tlik] coupling;
  matrix[N,Tlik] convout;

  // Ravg ~ normal(1.0,1.0);
  coupling_rate ~ normal(0.0, .1);
  dispersion ~ normal(0.0,10.0);

  // GP prior density
  eta ~ std_normal();
  gp_length_scale ~ normal(0.0,5.0);
  gp_sigma ~ normal(0.0, 0.5);
  local_sigma ~ normal(0.0, 0.5);
  global_sigma ~ normal(0.0, 1.0);


  // metapopulation infection rate model
  convout = METAPOP_metapop(Rt,convlik,coupling_rate);

  /* old coupling code
  for (i in 1:Tlik) {
    coupling[i] = 0.0;
    for (j in 1:N)
      coupling[i] += Rt[j] * convlik[j][i];
    coupling[i] = coupling[i] / N;
  }

  for (j in 1:N) {
    for (i in 1:Tlik) {
      Count[j,Tcond+i] ~ neg_binomial_2(
          (1.0-coupling_rate) *  Rt[j] * convlik[j,i] +
          coupling_rate * coupling[i], dispersion);
    }
  }
  */

  // compute likelihoods
  for (j in 1:N) 
    for (i in 1:Tlik) 
      Count[j,Tcond+i] ~ LIKELIHOOD_likelihood(convout[j,i], dispersion);

}

generated quantities {
  real Ravg = mean(Rt);
  vector[Tpred] Ppred[N];
  vector[Tproj] Cproj[N]; // faster than matrix

  // predictive probability of future counts
  for (i in 1:Tpred) {
    real convavg;
    convavg = 0.0;
    for (j in 1:N) 
      convavg += Rt[j] * convpred[j,i]; 
    convavg = convavg / N;

    for (j in 1:N) {
      Ppred[j,i] = exp(neg_binomial_2_lpmf(Count[j,Tcur+i] |
          (1.0-coupling_rate) *  Rt[j] * convpred[j,i] +
          coupling_rate * convavg, 
          dispersion));
    }
  }
  // forecasting *mean* counts given parameters
  for (i in 1:Tproj) {
    real convavg;
    vector[N] convprojall;
    convavg = 0.0;
    for (j in 1:N) {
      convprojall[j] = convproj[j,i] + 
          dot_product(Cproj[j][1:(i-1)], infprofile_rev[D-i+2:D]);
      convavg += Rt[j] * convprojall[j];
    }
    convavg = convavg / N;

    for (j in 1:N) {
      // Cproj[j,i] = neg_binomial_2_rng( 
      //     (1.0-coupling_rate) *  Rt[j] * convprojall[j] +
      //     coupling_rate * convavg, 
      //     dispersion);
      Cproj[j,i] = ( 
          (1.0-coupling_rate) *  Rt[j] * convprojall[j] +
          coupling_rate * convavg);
    }

  }

}


