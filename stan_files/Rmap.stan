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

  real local_var(real local_sigma2) {
    return local_sigma2;
  }

  real global_var(real global_sigma2) {
    return global_sigma2;
  }

  real none_var(real param) {
    return 0.000000001; // a small but positive number for numerical stability
  }

  // Meta-population infection rate model choices

  matrix uniform_out_metapop(
      vector Rt, matrix convlik, real coupling_rate, matrix flux, matrix convflux) {
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

  matrix uniform_in_metapop(
      vector Rt, matrix convlik, real coupling_rate, matrix flux, matrix convflux) {
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

  matrix radiation_out_metapop(
      vector Rt, matrix convlik, real coupling_rate, matrix flux, matrix convflux) {
    int T = cols(convlik);
    int N = rows(convlik);
    matrix[N,T] convout;
    
    for (i in 1:T) {
      for (j in 1:N) {
        convout[j,i] = Rt[j] * (1.0-coupling_rate) * convlik[j,i];
        for (k in 1:N) {
          convout[j,i] += Rt[k] * coupling_rate * convlik[k,i] * flux[k,j];
        } 
      }
    }
    return convout;
  }

  matrix radiation_in_metapop(
      vector Rt, matrix convlik, real coupling_rate, matrix flux, matrix convflux) {
    int T = cols(convlik);
    int N = rows(convlik);
    matrix[N,T] convout;
    
    for (i in 1:T) {
      for (j in 1:N) {
        convout[j,i] = Rt[j] * (
            (1.0-coupling_rate) * convlik[j,i] + 
            coupling_rate * convflux[j,i]
        );
      }
    }
    return convout;
  }


  matrix none_metapop(
      vector Rt, matrix convlik, real coupling_rate, matrix flux, matrix convflux) {
    int T = cols(convlik);
    int N = rows(convlik);
    matrix[N,T] convout;

    for (j in 1:N)
      for (i in 1:T)
        convout[j,i] = Rt[j] * convlik[j,i];

    return convout;
  }
  matrix compute_flux(matrix convlik, matrix flux) {
    int T = cols(convlik);
    int N = rows(convlik);
    matrix[N,T] convflux;
    for (i in 1:T) {
      for (j in 1:N) {
        convflux[j,i] = 0.0;
        for (k in 1:N) {
          convflux[j,i] += convlik[k,i] * flux[k,j];
        }
      }
    }
    return convflux;
  }

  // Case count likelihood choices

  real poisson_likelihood_lpmf(int count, real mu, real precision) {
    return poisson_lpmf(count | mu);
  }

  real negative_binomial_2_likelihood_lpmf(int count, real mu, real precision) {
    return neg_binomial_2_lpmf(count | mu, precision);
  }

  real negative_binomial_3_likelihood_lpmf(int count, real mu, real x) {
    return neg_binomial_2_lpmf(count | mu, mu / x);
  }


  real gig_lpdf(real x, int p, real a, real b) {
    return p * 0.5 * log(a / b)
      - log(2 * modified_bessel_second_kind(p, sqrt(a * b)))
      + (p - 1) * log(x)
      - (a * x + b / x) * 0.5;
  }
}

///////////////////////////////////////////////////////////////////////////

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
  matrix[N,N] flux;         // fluxes for radiation metapopulation model
}

///////////////////////////////////////////////////////////////////////////

transformed data {
  int Tcur = Tcond+Tlik;    // index of day on which we are estimating Rt
  int Tpred = Tall-Tcur;    // number of days to calculate predictive probabilities for

  vector[Tall] Creal[N];     // real type version of Count
  vector[D] infprofile_rev; // reversed infection profile

  // precompute convolutions between Count and infprofile 
  matrix[N,Tlik] convlik;      // for use in likelihood computation
  matrix[N,Tlik] convlikflux;  // for use in likelihood computation
  matrix[N,Tpred] convpred;    // for use in predictive probs of future counts
  matrix[N,Tpred] convpredflux;// for use in predictive probs of future counts
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
  convlikflux = compute_flux(convlik,flux);
  convpredflux = compute_flux(convpred,flux);
}

///////////////////////////////////////////////////////////////////////////

parameters {
  real<lower=0> gp_length_scale;
  real<lower=0> gp_sigma;
  real<lower=0> global_sigma;
  real<lower=0> local_sigma2[N];
  real<lower=0> local_scale;
  vector[N] eta;

  real<lower=0> precision;
  // real<lower=0> Ravg;
  real<lower=0,upper=1> coupling_rate;
}


///////////////////////////////////////////////////////////////////////////

transformed parameters {
  vector[N] Rt;                 // instantaneous reproduction number

  {
    matrix[N,N] K;
    matrix[N,N] L;
    real global_sigma2 = square(global_sigma);

    // GP prior.
    K = SPATIAL_kernel(geodist, gp_sigma, gp_length_scale); // kernel
    for (i in 1:N) {
      K[i,i] = K[i,i] + LOCAL_var(local_sigma2[i]);
    }
    K = K + GLOBAL_var(global_sigma2);

    L = cholesky_decompose(K);
    Rt = exp(L * eta);
  }
}


///////////////////////////////////////////////////////////////////////////

model {
  vector[Tlik] coupling;
  matrix[N,Tlik] convout;

  // Ravg ~ normal(1.0,1.0);
  coupling_rate ~ normal(0.0, .5);
  precision ~ normal(0.0,10.0);

  // GP prior density
  eta ~ std_normal();
  // gp_length_scale ~ normal(0.0,5.0);
  gp_length_scale ~ gig(2, 2.0, 2.0);
  gp_sigma ~ normal(0.0, 0.5);
  global_sigma ~ normal(0.0, 0.5);
  local_scale ~ normal(0.0, 0.1);
  for (i in 1:N) {
    local_sigma2[i] ~ exponential(0.5 / square(local_scale));
  }
  // metapopulation infection rate model
  convout = METAPOP_metapop(Rt,convlik,coupling_rate,flux,convlikflux);

  // compute likelihoods
  for (j in 1:N) 
    for (i in 1:Tlik) 
      Count[j,Tcond+i] ~ OBSERVATION_likelihood(convout[j,i], precision);

}



///////////////////////////////////////////////////////////////////////////

generated quantities {
  real R0;
  matrix[N,Tpred] Ppred;
  matrix[N,Tproj] Cproj; 

  // Estimated R0 over all areas
  {
    real convsum = 0.0;
    R0 = 0.0;
    for (i in 1:Tlik) {
      for (j in 1:N) {
        R0 += Rt[j] * convlik[j,i];
        convsum += convlik[j,i];
      }
    }
    R0 = R0 / convsum;
  }

  // predictive probability of future counts
  {
    matrix[N,Tpred] convout = METAPOP_metapop(Rt,convpred,coupling_rate,flux,convpredflux);
    for (i in 1:Tpred)
      for (j in 1:N)
        Ppred[j,i] = exp(OBSERVATION_likelihood_lpmf(Count[j,Tcur+i] |
            convout[j,i], precision));
  }

  // forecasting *mean* counts given parameters
  {
    matrix[N,1] convprojall;
    matrix[N,1] convprojflux;
    for (i in 1:Tproj) {
      for (j in 1:N) 
        convprojall[j,1] = convproj[j,i] + 
            dot_product(Cproj[j][1:(i-1)], infprofile_rev[(D-i+2):D]);
      convprojflux = compute_flux(convprojall,flux);
      Cproj[:,i] = METAPOP_metapop(Rt,convprojall,coupling_rate,flux,convprojflux)[:,1];
    }
  }

}


