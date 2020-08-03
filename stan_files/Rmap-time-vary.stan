functions {
  // GP kernel choices

  matrix exp_quad_kernel(matrix dist, real gp_sigma, real gp_length_scale) {
    return square(gp_sigma) * exp( (-1.0 / (2.0 * square(gp_length_scale))) * dist .* dist);
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

  real[,,] uniform1_metapop(
      matrix Rt, real[,,] convlik, real coupling_rate) {
    // Rt N*M
    // convlik [M,N,T]

    int M = size(convlik);
    int N = size(convlik[1]);
    int T = size(convlik[1,1]);

    real convavg[M,T];
    real convout[M,N,T];
    
    convavg = rep_array(0.0, M,T);

    // accessing in optimal order. See https://mc-stan.org/docs/2_23/stan-users-guide/indexing-efficiency-section.html
    for (k in 1:M) {
      for (j in 1:N) {
        for (i in 1:T) {
          convavg[k,i] += convlik[k,j,i] / N;
        }
      }
    }

    for (k in 1:M)
      for (j in 1:N)
        for (i in 1:T)
          convout[k,j,i] = (
              (1.0-coupling_rate) *  Rt[j, k] * convlik[k,j,i] +
              coupling_rate * convavg[k,i]
          );
    return convout;
  }

  real[,,] uniform2_metapop(
      matrix Rt, real[,,] convlik, real coupling_rate) {
    // Rt N*M
    // convlik [M,N,T]

    int M = size(convlik);
    int N = size(convlik[1]);
    int T = size(convlik[1,1]);

    real convavg[M,T];
    real convout[M,N,T];
    
    convavg = rep_array(0.0, M,T);

    // accessing in optimal order. See https://mc-stan.org/docs/2_23/stan-users-guide/indexing-efficiency-section.html
    for (k in 1:M) {
      for (j in 1:N)
        for (i in 1:T) {
          convavg[k,i] += convlik[k,j,i] / N;
      }
    }
    // convavg = convavg / N;

    for (k in 1:M)
      for (j in 1:N)
        for (i in 1:T)
          convout[k,j,i] = Rt[j,k] * (
              (1.0-coupling_rate) * convlik[k,j,i] +
              coupling_rate * convavg[k,i]
          );

    return convout;
  }

  real[,,] none_metapop(
      matrix Rt, real[,,] convlik, real coupling_rate) {
    // Rt N*M
    // convlik [M,N,T]

    int M = size(convlik);
    int N = size(convlik[1]);
    int T = size(convlik[1,1]);

    real convout[M,N,T];

    // accessing in optimal order. See https://mc-stan.org/docs/2_23/stan-users-guide/indexing-efficiency-section.html
    for (k in 1:M)
      for (j in 1:N)
        for (i in 1:T)
          convout[k,j,i] = Rt[j, k] * convlik[k,j,i];

    return convout;
  }

  // Case count likelihood choices

  real poisson_likelihood_lpmf(int count, real mean, real dispersion) {
    return poisson_lpmf(count | mean);
  }

  real negative_binomial_likelihood_lpmf(int count, real mean, real dispersion) {
    return neg_binomial_2_lpmf(count | mean, dispersion);
  }

  real gig_lpdf(real x, int p, real a, real b) {
    return p * 0.5 * log(a / b)
      - log(2 * modified_bessel_second_kind(p, sqrt(a * b)))
      + (p - 1) * log(x)
      - (a * x + b / x) * 0.5;
  }
}

data {
  int<lower=1> N;           // number of regions
  int<lower=1> M;           // number of time steps
  int<lower=1> D;           // length of infection profile
  int<lower=1> Tall;        // number of all days in case count time series
  int<lower=1> Tcond;       // number of days we will condition on
  int<lower=1> Tlik;        // number of days for likelihood computation
  int<lower=0> Tproj;       // number of days to forecast
  int<lower=0> Tstep;       // number of days to step for each time step of Rt prediction
  int Count[N, Tall];       // case counts,
  // vector[2] geoloc[N];   // geo locations of regions
  vector[D] infprofile;     // infection profile aka serial interval distribution
  matrix[N,N] geodist;      // distance between locations
  matrix[M,M] timedist;     // distance between time samples
}

transformed data {
  int Tcur = Tcond+Tlik;    // index of day on which we are estimating Rt
  int Tpred = Tall-Tcur;    // number of days to calculate predictive probabilities for

  vector[Tall] Creal[N];     // real type version of Count
  vector[D] infprofile_rev; // reversed infection profile

  // precompute convolutions between Count and infprofile 
  real convlik[M,N,Tlik];   // for use in likelihood computation. Batches of N places repeated for M time segments (i.e. [N=1 M=1, N=2 M=1, ..., N=N M=1, N=1 M=2, ...])
  real convpred[1,N,Tpred];    // for use in predictive probs of future counts. Extra dimension at front for cosistency (makes function definition much cleaner)
  real convproj[1,N,Tproj];    // for use in forecasting into future. Extra dimension at front for cosistency (makes function definition much cleaner)

  // reverse infection profile
  for (i in 1:D)
    infprofile_rev[i] = infprofile[D-i+1];

  for (j in 1:N) {
    for (i in 1:Tall)
      Creal[j,i] = Count[j,i];

    // precompute convolutions between counts and infprofile
    // compute for each time offset - NOTE: not the fastest ordering of these access options (see https://mc-stan.org/docs/2_23/stan-users-guide/indexing-efficiency-section.html), but assuming okay as this is only done once
    for (k in 1:M){
      for (i in 1:Tlik) {
        int L = min(D,Tcond+i-1-((M-k) * Tstep)); // length of infection profile that overlaps with case counts 
        convlik[k,j,i] = dot_product(Creal[j][Tcond+i-((M-k) * Tstep)-L:Tcond-1+i-((M-k) * Tstep)], infprofile_rev[D-L+1:D]);
      }
    }
    
    for (i in 1:Tpred) {
      int L = min(D,Tcur+i-1); // length of infection profile that overlaps with case counts 
      convpred[1,j,i] = dot_product(Creal[j][Tcur-L+i:Tcur-1+i], infprofile_rev[D-L+1:D]);
    }
    for (i in 1:Tproj) {
      int L = min(D,Tcur+i-1); // length of infection profile that overlaps with case counts 
      convproj[1,j,i] = dot_product(Creal[j][Tcur-L+i:Tcur], infprofile_rev[D-L+1:D-i+1]);
    }
  }
}

parameters {
  real<lower=0> gp_space_length_scale;
  real<lower=0> gp_space_sigma;

  real<lower=0> gp_time_length_scale;
  real<lower=0> gp_time_sigma;
  
  real<lower=0> local_sigma;
  real<lower=0> global_sigma;

  vector[N*M] eta;
  vector[N*M] epislon;

  real<lower=0> dispersion;
  // real<lower=0> Ravg;
  real<lower=0,upper=1> coupling_rate;
}

transformed parameters {
  matrix[N, M] Rt;                 // instantaneous reproduction number.

  {
    matrix[N,N] K_space;
    matrix[M,M] K_time;
    matrix[N,N] L_space;
    matrix[M,M] L_time;
    real local_sigma2 = square(local_sigma);
    real global_sigma2 = square(global_sigma);

    // GP prior.
    K_space = SPATIAL_kernel(geodist, gp_space_sigma, gp_space_length_scale); // space kernel
    K_time  = TEMPORAL_kernel(timedist, gp_time_sigma, gp_time_length_scale); // time kernel

    for (i in 1:N) {
      K_space[i,i] = K_space[i,i] + LOCAL_var(local_sigma2);
    }
    K_space = K_space + GLOBAL_var(global_sigma2);

    L_space = cholesky_decompose(K_space);
    L_time = cholesky_decompose(K_time);

    // Compute (K_time (*) K_space) * eta via efficient kronecker trick. Don't reshape back to vector for convinience.
    Rt = exp(L_space * to_matrix(eta, N, M) * L_time');
  }
}

model {
  vector[Tlik] coupling;
  real convout[M,N,Tlik];

  // Ravg ~ normal(1.0,1.0);
  coupling_rate ~ normal(0.0, .1);
  dispersion ~ normal(0.0,10.0);

  // GP prior density
  eta ~ std_normal();
  epislon ~ std_normal();

  gp_space_length_scale ~ gig(2, 2.0, 2.0);
  gp_space_sigma ~ normal(0.0, 0.5);

  gp_time_length_scale ~ gig(7, 1.0, 1.0);
  gp_time_sigma ~ normal(0.0, 0.5);

  local_sigma ~ normal(0.0, 0.5);
  global_sigma ~ normal(0.0, 1.0);

  // metapopulation infection rate model
  convout = METAPOP_metapop(Rt,convlik,coupling_rate);

  // compute likelihoods
  for (k in 1:M)
    for (i in 1:Tlik) 
      for (j in 1:N) 
        Count[j,Tcond+i-((M-k)*Tstep)] ~ OBSERVATION_likelihood(convout[k,j,i], dispersion);

}

generated quantities {
  row_vector[M] R0;
  real Ppred[N,Tpred];
  matrix[N,Tproj] Cproj; 

  // Estimated R0 over all areas
  {
    for (k in 1:M){
      real convsum = 0.0;
      R0[k] = 0.0;
      for (i in 1:Tlik) {
        for (j in 1:N) {
          R0[k] += Rt[j, k] * convlik[k,j,i];
          convsum += convlik[k,j,i];
        }
      }
      R0[k] = R0[k] / convsum;
    }
  }

  // predictive probability of future counts
  {
    real convout[1,N,Tpred] = METAPOP_metapop(Rt,convpred,coupling_rate);
    for (j in 1:N)
      for (i in 1:Tpred)
        Ppred[j,i] = exp(OBSERVATION_likelihood_lpmf(Count[j,Tcur+i] |
            convout[1,j,i], dispersion));
  }

  // forecasting *mean* counts given parameters
  {
    real convprojall[1,N,1];
    for (i in 1:Tproj) {
      for (j in 1:N) 
        convprojall[1,j,1] = convproj[1,j,i] + 
            dot_product(Cproj[j,1:(i-1)], infprofile_rev[D-i+2:D]);
      Cproj[:,i] = to_vector(METAPOP_metapop(to_matrix(Rt[:,M]),convprojall,coupling_rate)[1,:,1]);
    }
  }

}


