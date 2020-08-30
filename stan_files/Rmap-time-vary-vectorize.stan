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
    // return rep_matrix(0.0, rows(dist), cols(dist));
    return diag_matrix(rep_vector(0.000000001, cols(dist))); // very small diagonal to allow cholesky factorisation
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

  matrix[] none_metapop(matrix Rin, matrix[] convlik) {
    int N = rows(Rin);
    int M = cols(Rin);
    int T = cols(convlik[1]);

    matrix[N,T] convin[M];

    for (m in 1:M) {
      convin[m] = diag_pre_multiply(col(Rin, m), convlik[m]);
    }

    return convin;
  }

  matrix[] in_compute_flux(matrix[] convlik, matrix fluxt) {
    int M = size(convlik);
    int N = rows(convlik[1]);
    int T = cols(convlik[1]);
    int F1 = rows(fluxt);

    matrix[F1, N*T] flux_in[M];
    matrix[F1*N, N] reshaped_fluxt;
    reshaped_fluxt = to_matrix(fluxt, F1*N, N);

    for (m in 1:M){
      flux_in[m] = to_matrix(reshaped_fluxt * convlik[m], F1, N*T);
    }

    return flux_in;
  }

  matrix[] in_metapop(matrix Rin, matrix[] convlik, 
      row_vector[] fluxproportions, matrix[] convflux) {
    // int M = cols(Rin);
    int M = size(convlik);
    int N = rows(convlik[1]);
    int T = cols(convlik[1]);

    matrix[N,T] convin[M];

    for (m in 1:M) {
      convin[m] = diag_pre_multiply(col(Rin, m), to_matrix(fluxproportions[m] * convflux[m], N, T));
    }

    return convin;
  }

  matrix[] in_out_metapop(
      matrix Rin, matrix Rout, matrix[] convlik, row_vector[] fluxproportions, matrix fluxt) {
    int M = size(convlik);
    int N = rows(convlik[1]);
    int T = cols(convlik[1]);
    int F = rows(fluxt);
    // convin =                   [M][Nin,T] 
    // diag_pre_multiply(         [M][Nin,T]
    //   Rin,                     [M][Nin]
    //   (
    //     to_matrix(
    //       (                
    //         fluxproportions    [M][F] 
    //         *
    //         fluxt              [F,Nin*Nout] 
    //       )                    [M][Nin*Nout]
    //       , Nin, Nout
    //     )                      [M][Nin,Nout]
    //     *
    //     diag_pre_multiply(
    //       Rout,                [M][Nout] 
    //       convlik              [M][Nout,T]
    //     )                      [M][Nout,T]
    //   )                        [M][Nin,T]
    // );                         [M][Nin,T]
    matrix[N,T] convin[M];
    for (m in 1:M) {
      convin[m] = diag_pre_multiply(col(Rin, m), 
        to_matrix(fluxproportions[m] * fluxt, N, N) *
        diag_pre_multiply(col(Rout, m),convlik[m]));
    }
    return convin;
  }

  matrix[] metapop(
      int do_metapop, int do_in_out,
      matrix Rin, matrix Rout, matrix[] convlik, matrix[] convflux,
      row_vector[] fluxproportions, matrix fluxt) {
    if (do_metapop) {
      if (do_in_out) {
        return in_out_metapop(Rin,Rout,convlik,fluxproportions,fluxt);
      } else {
        return in_metapop(Rin,convlik,fluxproportions,convflux);
      }
    } else {
      return none_metapop(Rin,convlik);
    }
  }
  

  // Case count likelihood choices

  real poisson_likelihood_lpmf(int count, real mean, real precision) {
    return poisson_lpmf(count | mean);
  }

  real negative_binomial_2_likelihood_lpmf(int count, real mu, real precision) {
    return neg_binomial_2_lpmf(count | mu, precision);
  }

  real negative_binomial_3_likelihood_lpmf(int count, real mu, real x) {
    return neg_binomial_2_lpmf(count | mu, mu / x);
  }

  // generalized inverse gamma distribution

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
  matrix[M,M] timecorcut;   // matrix specifying which time points should be correlated (to account for lockdown)
  int<lower=0,upper=1> do_metapop;
  int<lower=0,upper=1> do_in_out;
  int F;
  matrix[N,N] flux[F];      // fluxes for radiation metapopulation model
}

transformed data {
  int Tcur = Tcond+Tlik;    // index of day on which we are estimating Rt
  int Tpred = Tall-Tcur;    // number of days to calculate predictive probabilities for
  int F1 = F+1;
  vector[max(1,F)] ones = rep_vector(1.0,max(1,F));

  vector[Tall] Creal[N];     // real type version of Count
  vector[D] infprofile_rev; // reversed infection profile

  // precompute convolutions between Count and infprofile 
  matrix[N,Tlik] convlik[M];      // for use in likelihood computation
  matrix[N,Tpred] convpred[1];    // for use in predictive probs of future counts
  matrix[N,Tproj] convproj[1];    // for use in forecasting into future 

  matrix[F1,N*N] fluxt;      // transposed flux matrices
  matrix[F1,N*Tlik] convlikflux[M];
  matrix[F1,N*Tpred] convpredflux[1];

  // reverse infection profile
  for (i in 1:D)
    infprofile_rev[i] = infprofile[D-i+1];

  {
    fluxt[1,] = to_row_vector(diag_matrix(rep_vector(1.0,N)));
    for (f in 2:F1)
      fluxt[f,] = to_row_vector(flux[f-1]');
  }

  for (j in 1:N)
    for (i in 1:Tall)
      Creal[j,i] = Count[j,i];

  // precompute convolutions between counts and infprofile
  // compute for each time offset - NOTE: not the fastest ordering of these access options (see https://mc-stan.org/docs/2_23/stan-users-guide/indexing-efficiency-section.html), but assuming okay as this is only done once
  for (j in 1:N) {
    for (k in 1:M) {
      for (i in 1:Tlik) {
        int L = min(D,Tcond+i-1-((M-k) * Tstep)); // length of infection profile that overlaps with case counts 
        convlik[k,j,i] = dot_product(Creal[j][Tcond+i-((M-k) * Tstep)-L:Tcond-1+i-((M-k) * Tstep)], infprofile_rev[D-L+1:D])+1e-6;
      }
    }
    
    for (i in 1:Tpred) {
      int L = min(D,Tcur+i-1); // length of infection profile that overlaps with case counts 
      convpred[1,j,i] = dot_product(Creal[j][Tcur-L+i:Tcur-1+i], infprofile_rev[D-L+1:D])+1e-6;
    }
    for (i in 1:Tproj) {
      int L = min(D,Tcur+i-1); // length of infection profile that overlaps with case counts 
      convproj[1,j,i] = dot_product(Creal[j][Tcur-L+i:Tcur], infprofile_rev[D-L+1:D-i+1])+1e-6;
    }
  }

  if (do_metapop && !do_in_out) {
    convlikflux = in_compute_flux(convlik,fluxt);
    convpredflux = in_compute_flux(convpred,fluxt);
  }
}

parameters {
  real<lower=0> gp_space_length_scale;
  real<lower=0> gp_space_sigma;

  real<lower=0> gp_time_length_scale;
    
  matrix<lower=0>[N,M] local_exp;
  real<lower=0> local_scale;
  real<lower=0> global_sigma;

  vector[N*M] eta_in;
  vector[N*M] eta_out;

  vector[N*M] epsilon_in;
  vector[N*M] epsilon_out;

  real<lower=0> precision;
  // real<lower=0> Ravg;
  real<lower=0,upper=1> coupling_rate[M];
  simplex[max(1,F)] flux_probs;
}

transformed parameters {
  matrix[N, M] Rin;                 // instantaneous reproduction number
  matrix[N, M] Rout;                 // instantaneous reproduction number
  matrix[N, M] local_sigma;
  row_vector[F1] fluxproportions[M];
  matrix[N,Tlik] convlikout[M];
  {
    matrix[N,N] K_space;
    matrix[M,M] K_time;
    matrix[N,N] L_space;
    matrix[M,M] L_time;
    real global_sigma2 = square(global_sigma);

    // GP space kernel
    K_space = SPATIAL_kernel(geodist, gp_space_sigma, gp_space_length_scale); // space kernel
    K_space = K_space + GLOBAL_var(global_sigma2); // Add global noise

    // GP time kernel
    K_time  = TEMPORAL_kernel(timedist, 1.0, gp_time_length_scale);
    K_time  = K_time .* timecorcut;  // Zero out uncorrelated time entries

    L_space = cholesky_decompose(K_space);
    L_time = cholesky_decompose(K_time);

    // Noise kernel (reshaped from (N*M, N*M) to (N,M) since diagonal)
    for (m in 1:M) {
      for (n in 1:N) {
        local_sigma[n, m] = sqrt(LOCAL_var(local_exp[n, m] * (2 * square(local_scale))));
      }
    }
    
    // Compute (K_time (*) K_space) * eta via efficient kronecker trick. Don't reshape back to vector for convinience.
    // Add on the location-time dependant noise as well
    Rin = exp((L_space * to_matrix(eta_in, N, M) * L_time') + (local_sigma .* to_matrix(epsilon_in, N, M)));
    if (do_metapop && do_in_out) {
      Rout = exp((L_space * to_matrix(eta_out, N, M) * L_time') + (local_sigma .* to_matrix(epsilon_out, N, M)));
    } else {
      Rout = rep_matrix(1.0,N,M);
    }

    // metapopulation infection rate model
    if (do_metapop) {
      for (m in 1:M){
        fluxproportions[m, 1] = 1.0-coupling_rate[m];
      for (f in 2:F1)
        fluxproportions[m, f] = coupling_rate[m]*(1.0-flux_probs[f-1]);
      }
    } else {
      for (m in 1:M){
        fluxproportions[m, 1] = 1.0;
        for (f in 2:F1)
          fluxproportions[m, f] = 0.0;
      }
      
    }
    convlikout = metapop(do_metapop,do_in_out,
        Rin,Rout,convlik,convlikflux,fluxproportions,fluxt);
  }
}

model {
  coupling_rate ~ normal(0.0, .25);
  flux_probs ~ dirichlet(ones);
  precision ~ normal(0.0,5.0);

  // GP prior density
  eta_in ~ std_normal();
  eta_out ~ std_normal();
  epsilon_in ~ std_normal();
  epsilon_out ~ std_normal();

  gp_space_length_scale ~ gig(5, 5.0, 5.0);
  gp_space_sigma ~ normal(0.0, 0.25);

  gp_time_length_scale ~ gig(14, 1.0, 1.0);

  local_scale ~ normal(0.0, 0.5);
  for (j in 1:M){
    for (i in 1:N) {
      local_exp[i, j] ~ exponential(1.0);
      // local_sigma2[i, j] ~ exponential(0.5 / square(local_scale)); // reparameterised
    }
  }
  global_sigma ~  normal(0.0, 0.25);

  // compute likelihoods
  for (k in 1:M) {
    for (i in 1:Tlik) {
      for (j in 1:N) {
        Count[j,Tcond+i-((M-k)*Tstep)] ~ OBSERVATION_likelihood(convlikout[k,j,i], precision);
      }
    }
  }
}

generated quantities {
  real R0[M];
  vector[N] Rt[M];
  matrix[N,Tpred] Ppred;
  matrix[N,M*Tlik] Cpred;
  matrix[N,Tproj] Cproj; 

  // Estimated R0 and Rt for all areas
  {
    matrix[N, M] oneN = rep_matrix(1.0,N,M);
    vector[Tlik] oneT = rep_vector(1.0,Tlik);
    matrix[N,Tlik] convone[M] = metapop(do_metapop,do_in_out,
        oneN,oneN,convlik,convlikflux,fluxproportions,fluxt);
    for (m in 1:M) {
      R0[m] = sum(convlikout[m]) / sum(convone[m]);
      Rt[m] = (convlikout[m] * oneT) ./ (convone[m] * oneT);
    }
  }
  
  
  // predictive probability of future counts
  {
    row_vector[F1] pred_fluxproportions[1];
    matrix[N,1] pred_Rin;
    matrix[N,1] pred_Rout;
    matrix[N,Tpred] convpredout[1];
    
    pred_fluxproportions[1] = fluxproportions[M];
    pred_Rin = block(Rin, 1, M, N, 1);
    pred_Rout = block(Rout, 1, M, N, 1);

    convpredout = metapop(do_metapop,do_in_out,
        pred_Rin,pred_Rout,convpred,convpredflux,pred_fluxproportions,fluxt);

    for (i in 1:Tpred)
      for (j in 1:N)
        Ppred[j,i] = exp(OBSERVATION_likelihood_lpmf(Count[j,Tcur+i] |
            convpredout[1,j,i], precision));
  }


  // posterior predictive expected counts
  {
    matrix[N,Tlik] convlikout[M] = metapop(do_metapop,do_in_out,
        Rin,Rout,convlik,convlikflux,fluxproportions,fluxt);
    for (k in 1:M)
      Cpred[,(1+(k-1)*Tlik):(k*Tlik)] = convlikout[k];
  }


  // forecasting expected counts given parameters
  {
    matrix[N,1] convprojall[1];
    matrix[F1,N*1] convprojflux[1];
    row_vector[F1] proj_fluxproportions[1];
    matrix[N,1] proj_Rin;
    matrix[N,1] proj_Rout;

    proj_fluxproportions[1] = fluxproportions[M];
    proj_Rin = block(Rin, 1, M, N, 1);
    proj_Rout = block(Rout, 1, M, N, 1);

    for (i in 1:Tproj) {
      for (j in 1:N) 
        convprojall[1,j,1] = convproj[1,j,i] + 
            dot_product(Cproj[j,1:(i-1)], infprofile_rev[(D-i+2):D]);
      if (do_metapop && !do_in_out)
        convprojflux = in_compute_flux(convprojall,fluxt);
      Cproj[:,i] = metapop(do_metapop,do_in_out,
          proj_Rin,proj_Rout,convprojall,convprojflux,fluxproportions,fluxt)[1,:,1];
    }
  }
}

