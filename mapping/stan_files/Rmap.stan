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

  matrix kernel(int KERNEL, matrix dist, real gp_sigma, real gp_length_scale,
      int NONE_KERNEL, int EXP_QUAD_KERNEL,  
      int MATERN12_KERNEL, int MATERN32_KERNEL, int MATERN52_KERNEL) {
    if (KERNEL == EXP_QUAD_KERNEL) {
      return exp_quad_kernel(dist, gp_sigma, gp_length_scale);
    } else if (KERNEL == MATERN12_KERNEL) {
      return matern12_kernel(dist, gp_sigma, gp_length_scale);
    } else if (KERNEL == MATERN32_KERNEL) {
      return matern32_kernel(dist, gp_sigma, gp_length_scale);
    } else if (KERNEL == MATERN52_KERNEL) {
      return matern52_kernel(dist, gp_sigma, gp_length_scale);
    } else if (KERNEL == NONE_KERNEL) {
      return none_kernel(dist, gp_sigma, gp_length_scale);
    }
    reject("Unknown kernel");
    return none_kernel(dist, gp_sigma, gp_length_scale);
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
      int DO_METAPOP, int DO_IN_OUT,
      matrix Rin, matrix Rout, matrix[] convlik, matrix[] convflux,
      row_vector[] fluxproportions, matrix fluxt) {
    if (DO_METAPOP) {
      if (DO_IN_OUT) {
        return in_out_metapop(Rin,Rout,convlik,fluxproportions,fluxt);
      } else {
        return in_metapop(Rin,convlik,fluxproportions,convflux);
      }
    } else {
      return none_metapop(Rin,convlik);
    }
  }
  

  // Case count likelihood choices

  real poisson_likelihood_lpmf(int count, real mean, real dispersion) {
    return poisson_lpmf(count | mean);
  }

  real negative_binomial_2_likelihood_lpmf(int count, real mu, real dispersion) {
    return neg_binomial_2_lpmf(count | mu, 1.0 / dispersion);
  }

  real negative_binomial_3_likelihood_lpmf(int count, real mu, real dispersion) {
    return neg_binomial_2_lpmf(count | mu, mu / dispersion);
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
  int<lower=1> Mstep;           // number of time steps
  int<lower=1> Tall;        // number of all days in case count time series
  int<lower=1> Tcond;       // number of days we will condition on
  int<lower=0> Tstep;       // number of days to step for each time step of Rt prediction
  int<lower=0> Mproj;           // number of time steps to forecast

  int Count[N,Tall];           // case counts
  int Clean_recon[N,Tall];     // case counts
  matrix[N,Tall] Clean_latent; // cleaned case counts

  // vector[2] geoloc[N];   // geo locations of regions
  int<lower=1> Tip;         // length of infection profile
  vector[Tip] infprofile;   // infection profile aka serial interval distribution
  int<lower=1> Tdp;         // length of infection profile
  vector[Tdp] delayprofile; // infection profile aka serial interval distribution
  int F;
  matrix[N,N] flux[F];      // fluxes for radiation metapopulation model


  matrix[N,N] geodist;      // distance between locations
  matrix[Mstep+Mproj,Mstep+Mproj] timedist;     // distance between time samples
  matrix[Mstep+Mproj,Mstep+Mproj] timecorcut;   // matrix specifying which time points should be correlated (to account for lockdown)

  real<lower=0> gp_space_scale; // 50km
  real<lower=0> gp_time_scale; // 3 weeks
  real<lower=0> gp_space_decay_scale;
  real<lower=0> gp_time_decay_scale;
  real fixed_gp_space_length_scale; // If positive, use this, otherwise use prior
  real fixed_gp_time_length_scale; // If positive, use this, otherwise use prior
  real<lower=0> coupling_mu_scale;
  real<lower=0> coupling_sigma_scale;
  real<lower=0> coupling_alpha_scale;
  int<lower=1,upper=5> SPATIAL_KERNEL;
  int<lower=1,upper=5> TEMPORAL_KERNEL;
  int<lower=0,upper=1> LOCAL_KERNEL;
  int<lower=0,upper=1> GLOBAL_KERNEL;
  int<lower=0,upper=1> DO_METAPOP;
  int<lower=0,upper=1> DO_IN_OUT;
  int<lower=1,upper=4> OBSERVATION_DATA;
  int<lower=1,upper=4> OBSERVATION_MODEL;
}

transformed data {
  int Tlik = Mstep*Tstep;     // Number of days to use for the likelihood computation
  int Tproj = Mproj*Tstep;    // The number of days to project forward
  int Tcur = Tcond+Tlik;    // index of day on which we are estimating Rt
  int Tpred = Tall - Tcur;    // number of days to calculate predictive probabilities for
  int Mpred = (Tall - Tcur + Tstep - 1) / Tstep; // Enough Tsteps to cover whole Tpred period. Funky math to deal with Stans deficiencies
  int F1 = F+1;
  vector[max(1,F)] ones = rep_vector(1.0,max(1,F));

  matrix[N,Tall] Creal;     // real type version of Clean
  vector[Tip] infprofile_rev; // reversed infection profile
  vector[Tdp] delayprofile_rev; // reversed infection profile
  vector[Tdp] delayprofile_cum; // reversed infection profile

  // precompute convolutions between Counts and infprofile 
  matrix[N,Tstep] convlik[Mstep];      // for use in likelihood computation
  int Tforw = max(Tpred,Tproj);
  int Mforw = max(Mpred,Mproj);
  matrix[N,Tforw] convforw[1];    // for use in forecasting into future 

  matrix[F1,N*N] fluxt;      // transposed flux matrices
  matrix[F1,N*Tstep] convlikflux[Mstep];

  int NONE_KERNEL = 1;
  int EXP_QUAD_KERNEL = 2;
  int MATERN12_KERNEL = 3;
  int MATERN32_KERNEL = 4;
  int MATERN52_KERNEL = 5;

  // OBSERVATION_MODEL
  int POISSON = 1;
  int NEG_BINOMIAL_2 = 2;
  int NEG_BINOMIAL_3 = 3;
  int GAUSSIAN = 4;

  // OBSERVATION_DATA
  int COUNTS = 1;
  int CLEANED_LATENT = 2;
  int CLEANED_RECON = 3;
  int INFECTION_REPORTS = 4;

  { // infection and delay profiles
    real s = 0.0;
    for (i in 1:Tip)
      infprofile_rev[i] = infprofile[Tip-i+1];
    for (i in 1:Tdp) {
      s += delayprofile[i];
      delayprofile_rev[Tdp-i+1] = delayprofile[i];
      delayprofile_cum[i] = s;
    }
  }
  { // flux matrices
    fluxt[1,] = to_row_vector(diag_matrix(rep_vector(1.0,N)));
    for (f in 2:F1)
      fluxt[f,] = to_row_vector(flux[f-1]');
  }

  if (OBSERVATION_DATA==CLEANED_LATENT || OBSERVATION_DATA==INFECTION_REPORTS) {
    Creal = Clean_latent;
  } else if (OBSERVATION_DATA==CLEANED_RECON) {
    for (j in 1:N) 
      for (i in 1:Tall)
        Creal[j,i] = Clean_recon[j,i];
  } else {
    for (j in 1:N) 
      for (i in 1:Tall)
        Creal[j,i] = Count[j,i];
  }

  // precompute convolutions between counts and infprofile
  // compute for each time offset - NOTE: not the fastest ordering of these access options (see https://mc-stan.org/docs/2_23/stan-users-guide/indexing-efficiency-section.html), but assuming okay as this is only done once
  for (j in 1:N) {
    for (k in 1:Mstep) {
      for (i in 1:Tstep) {
        int L = min(Tip,Tcond+i-1+((k-1) * Tstep)); // length of infection profile that overlaps with case counts 
        convlik[k,j,i] = dot_product(Creal[j,Tcond+i+((k-1)*Tstep)-L:Tcond+i+((k-1)*Tstep)-1], infprofile_rev[Tip-L+1:Tip]) + 1e-6;
      }
    }
    
    for (i in 1:Tforw) {
      int L = min(Tip,Tcur+i-1); // length of infection profile that overlaps with case counts 
      convforw[1,j,i] = dot_product(Creal[j,Tcur-L+i:Tcur], infprofile_rev[Tip-L+1:Tip-i+1])+1e-6;
    }
  }

  if (DO_METAPOP && !DO_IN_OUT) {
    convlikflux = in_compute_flux(convlik,fluxt);
  }

  // Mpred = 0;
  // Tpred = 0;
  // Mforw = 0;
  // Tforw = 0;
  // Tproj = 0;
}

parameters {
  real<lower=0.0> gp_sigma;
  real<lower=0.0,upper=0.632> gp_space_decay;
  real<lower=0.0,upper=0.632> gp_time_decay;
    
  matrix<lower=0.0>[N,Mstep+Mforw] local_exp;
  real<lower=0.0> local_scale;
  real<lower=0.0> global_sigma;

  vector[N*(Mstep+Mforw)] eta_in;
  vector[N*(Mstep+Mforw)] eta_out;

  vector[N*(Mstep+Mforw)] epsilon_in;
  vector[N*(Mstep+Mforw)] epsilon_out;

  real<lower=0.0> dispersion;
  // real<lower=0> Ravg;

  real coupling_mu; // prior mean for sigmoid^-1(coupling_rate) 
  real<lower=0> coupling_sigma; // prior scale for sigmoid^-1(coupling_rate) 
  real<lower=0,upper=1> coupling_alpha1; // 1-autocorrelation for sigmoid^-1(coupling_rate) 
  vector[Mstep+Mforw] coupling_eta; // noise for sigmoid^-1(coupling_rate) 

  simplex[max(1,F)] flux_probs;
}

transformed parameters {
  matrix[N, Mstep+Mforw] Rin;                 // instantaneous reproduction number
  matrix[N, Mstep+Mforw] Rout;                 // instantaneous reproduction number
  matrix[N, Mstep+Mforw] local_sigma;
  real coupling_alpha = 1-coupling_alpha1; // autocorrelation for sigmoid^-1(coupling_rate) 
  real coupling_rate[Mstep+Mforw];
  row_vector[F1] fluxproportions[Mstep+Mforw];
  matrix[N,Tstep] convlikout[Mstep];
  real gp_space_length_scale;
  real gp_time_length_scale;
  { // AR1 process for sigmoid^-1(coupling_rate)
    vector[Mstep+Mforw] Zt = rep_vector(0.0,Mstep+Mforw);
    real delta = sqrt(1.0-coupling_alpha);
    Zt[1] += coupling_eta[1];
    for (i in 2:(Mstep+Mforw))
      Zt[i] += coupling_alpha * Zt[i-1] + delta * coupling_eta[i];
    coupling_rate = to_array_1d(inv_logit(coupling_mu + coupling_sigma * Zt));
  }
  {
    matrix[N,N] K_space;
    matrix[Mstep+Mforw,Mstep+Mforw] K_time;
    matrix[N,N] L_space;
    matrix[Mstep+Mforw,Mstep+Mforw] L_time;
    real global_sigma2 = square(global_sigma);

    if (fixed_gp_space_length_scale<=0.0) {
      gp_space_length_scale = - gp_space_scale / log(1.0-gp_space_decay);
    } else {
      gp_space_length_scale = fixed_gp_space_length_scale;
    }
    if (fixed_gp_time_length_scale<=0.0) {
      gp_time_length_scale = - gp_time_scale / log(1.0-gp_time_decay);
    } else {
      gp_time_length_scale = fixed_gp_time_length_scale;
    }
 
    // GP space kernel
    K_space = kernel(SPATIAL_KERNEL,geodist, gp_sigma, gp_space_length_scale,
        NONE_KERNEL,EXP_QUAD_KERNEL,MATERN12_KERNEL,MATERN32_KERNEL,MATERN52_KERNEL);
    if (GLOBAL_KERNEL) 
      K_space = K_space + global_sigma2; // Add global noise

//print("K_space");
//print(K_space);

    // GP time kernel

    K_time  = kernel(TEMPORAL_KERNEL,timedist, 1.0, gp_time_length_scale,
        NONE_KERNEL,EXP_QUAD_KERNEL,MATERN12_KERNEL,MATERN32_KERNEL,MATERN52_KERNEL);
    K_time  = K_time .* timecorcut;  // Zero out uncorrelated time entries

//print("K_time");
//print(K_time);

    L_space = cholesky_decompose(K_space);
    L_time = cholesky_decompose(K_time);

    // Noise kernel (reshaped from (N*Mstep, N*Mstep) to (N,Mstep) since diagonal)
    if (LOCAL_KERNEL) {
      for (m in 1:(Mstep+Mforw)) {
        for (n in 1:N) {
          local_sigma[n, m] = sqrt((2 * square(local_scale)) * local_exp[n, m]);
        }
      }
    } else {
      local_sigma = rep_matrix(0.0,N,(Mstep+Mforw));
    }
    
//print("local_sigma");
//print(local_sigma);

    // Compute (K_time (*) K_space) * eta via efficient kronecker trick. 
    // Don't reshape back to vector for convinience.
    // Add on the location-time dependant noise as well
    Rin = exp(
        (L_space * to_matrix(eta_in, N, Mstep+Mproj) * L_time') + 
        (local_sigma .* to_matrix(epsilon_in, N, Mstep+Mproj))
    );
    if (DO_METAPOP && DO_IN_OUT) {
      Rout = exp(
          (L_space * to_matrix(eta_out, N, Mstep+Mproj) * L_time') + 
          (local_sigma .* to_matrix(epsilon_out, N, Mstep+Mproj))
      );
    } else {
      Rout = rep_matrix(1.0,N,Mstep+Mproj);
    }

    // metapopulation infection rate model
    if (DO_METAPOP) {
      for (m in 1:(Mstep+Mforw)){
        fluxproportions[m, 1] = 1.0-coupling_rate[m];
        for (f in 2:F1)
          fluxproportions[m, f] = coupling_rate[m]*flux_probs[f-1];
      }
    } else {
      for (m in 1:(Mstep+Mforw)){
        fluxproportions[m, 1] = 1.0;
        for (f in 2:F1)
          fluxproportions[m, f] = 0.0;
      }
      
    }
    convlikout= metapop(DO_METAPOP,DO_IN_OUT,
        block(Rin, 1, 1, N, Mstep),block(Rout, 1, 1, N, Mstep),convlik,convlikflux,fluxproportions[1:Mstep],fluxt);
  }
}

model {
  flux_probs ~ dirichlet(ones);
  dispersion ~ normal(0.0,5.0);

  // GP prior density
  eta_in ~ std_normal();
  eta_out ~ std_normal();
  epsilon_in ~ std_normal();
  epsilon_out ~ std_normal();


  gp_space_decay ~ normal(0.0,gp_space_decay_scale);
  gp_time_decay ~ normal(0.0,gp_time_decay_scale);
  gp_sigma ~ normal(0.0, 0.25);
  global_sigma ~ normal(0.0, 0.25);
  local_scale ~ normal(0.0, 0.1);
  for (j in 1:(Mstep+Mforw)){
    for (i in 1:N) {
      local_exp[i, j] ~ exponential(1.0);
      // local_sigma2[i, j] ~ exponential(0.5 / square(local_scale)); // reparameterised
    }
  }

  
  coupling_mu ~ normal(0.0, coupling_mu_scale);
  coupling_sigma ~ normal(0.0, coupling_sigma_scale);
  coupling_alpha1 ~ normal(0.0, coupling_alpha_scale);
  coupling_eta ~ std_normal();

  // compute likelihoods
  if (OBSERVATION_DATA == INFECTION_REPORTS) {
    matrix[N,Tcur] Einfect;
    Einfect[,1:Tcond] = Clean_latent[,1:Tcond];
    for (k in 1:Mstep)
      Einfect[,(Tcond+1+(k-1)*Tstep):(Tcond+k*Tstep)] = convlikout[k];
    for (s in 1:Tlik) {
      int t = Tcond+s;
      vector[N] Ereport = Einfect[,t-Tdp+1:t] * delayprofile_rev;
      for (j in 1:N)
        Count[j,t] ~ neg_binomial_2(
              Ereport[j],
              Ereport[j] / dispersion
        ); 
    }
  } else {
    for (k in 1:Mstep) {
      for (j in 1:N) {
        if (OBSERVATION_DATA == COUNTS) {

          reject("Currently invalid. Needs refactor")
          // if (OBSERVATION_MODEL == POISSON) {
          //   Count_lik_reduced[k,j] ~ poisson(
          //       convlikout_reduced[k,j,1] 
          //   );
          // } else if (OBSERVATION_MODEL == NEG_BINOMIAL_2) {
          //   Count_lik_reduced[k,j] ~ neg_binomial_2(
          //       convlikout_reduced[k,j,1],
          //       1.0 / dispersion
          //   );
          // } else if (OBSERVATION_MODEL == NEG_BINOMIAL_3) {
          //   Count_lik_reduced[k,j] ~ neg_binomial_2(
          //       convlikout_reduced[k,j,1],
          //       convlikout_reduced[k,j,1] / dispersion
          //   );
          // } else {
          //   reject("Invalid combination of OBSERVATION_DATA, OBSERVATION_MODEL found: ", OBSERVATION_DATA, ", ", OBSERVATION_MODEL);
          // }

        } else if (OBSERVATION_DATA == CLEANED_LATENT) {

          if (OBSERVATION_MODEL == GAUSSIAN) {
            for (i in 1:Tstep) {
              vector[2] loglik;
              real latent = Clean_latent[j,Tcond+i+(k-1)*Tstep];
              real Elatent = convlikout[k,j,i];
              loglik[1] = normal_lpdf( latent | Elatent, sqrt((1.0+dispersion)*Elatent) );
              loglik[2] = normal_lpdf( -latent | Elatent, sqrt((1.0+dispersion)*Elatent) );
              target += log_sum_exp(loglik);
            }
            // Clean_latent_lik_reduced[k,j] ~ normal(
            //     convlikout_reduced[k,j,1], 
            //     sqrt((1.0+dispersion)*convlikout_reduced[k,j,1])
            // );
          } else {
            reject("Invalid combination of OBSERVATION_DATA, OBSERVATION_MODEL found: ", OBSERVATION_DATA, ", ", OBSERVATION_MODEL);
          }

        } else if (OBSERVATION_DATA == CLEANED_RECON) {

          reject("Currently invalid. Needs refactor")
          // if (OBSERVATION_MODEL == POISSON) {
          //   Clean_recon_lik_reduced[k,j] ~ poisson(
          //       convlikout_reduced[k,j,1] 
          //   );
          // } else if (OBSERVATION_MODEL == NEG_BINOMIAL_2) {
          //   Clean_recon_lik_reduced[k,j] ~ neg_binomial_2(
          //       convlikout_reduced[k,j,1],
          //       1.0 / dispersion
          //   );
          // } else if (OBSERVATION_MODEL == NEG_BINOMIAL_3) {
          //   Clean_recon_lik_reduced[k,j] ~ neg_binomial_2(
          //       convlikout_reduced[k,j,1],
          //       convlikout_reduced[k,j,1] / dispersion
          //   );
          // } else {
          //   reject("Invalid combination of OBSERVATION_DATA, OBSERVATION_MODEL found: ", OBSERVATION_DATA, ", ", OBSERVATION_MODEL);
          // }

        }
      }
    }
  }
}

generated quantities {
  real Rt_all[Mstep+Mproj];
  vector[N] Rt[Mstep+Mproj];
  matrix[N,Tpred] Ppred;
  matrix[N,Mstep*Tstep] Cpred;
  matrix[N,Mproj*Tstep] Cproj;

  {
    matrix[N,Tcur+Tforw] Clatent;     // Latent epidemic values
    matrix[N,Tstep] convlik_forw[Mstep+Mforw]; // Infection pressure
    matrix[N,1] convlik_forw_reduced[Mstep+Mforw]; // for use in Rt computation
    matrix[N,Tstep] convlikout_forw[Mstep+Mforw]; // Infection pressure
    matrix[N,1] convlikout_forw_reduced[Mstep+Mforw]; // for use in Rt computation

    // Fill Clatent with observations up until the current time
    Clatent[,1:Tcur] = Creal[,1:Tcur];
    for (k in 1:Mstep) {
      convlikout_forw[k] = convlikout[k];
      convlik_forw[k] = convlik[k];
    }

    { // print stats
      if (uniform_rng(0.0,1.0)<1.0) {
        print(
          "space ", gp_space_length_scale,
          "; time ", gp_time_length_scale,
          "; sigmas gp ", gp_sigma," g ", global_sigma," l ", local_scale,
          "; dispersion ",dispersion
        );
      }
    }   

    // Forward simulate model and compute predictive probabilities
    if (Mforw > 0) {
      // Rollout stochasitc prediction of epidemic
      {
        matrix[F1,N*1] convforwflux[1];
        matrix[N,Mforw] forw_Rin;
        matrix[N,Mforw] forw_Rout;        

        // Pick out the posterior samples of the future parameters conditioned on past observations
        forw_Rin = block(Rin, 1, Mstep+1, N, Mforw);
        forw_Rout = block(Rout, 1, Mstep+1, N, Mforw);

        // For each timestep simulate the model forwards
        for (m in 1:Mforw) {
          for (t in 1:Tstep) {
            int i = (m-1) * Tstep + t;
            // Compute the infection pressure by adding on pressure from simulated counts
            for (j in 1:N) 
              convlik_forw[Mstep+m, j, t] = convforw[1,j,i] + 
                  dot_product(Clatent[j,Tcur+1:Tcur+(i-1)], infprofile_rev[(Tip-i+2):Tip]);

            // Compute flux
            if (DO_METAPOP && !DO_IN_OUT)
              convforwflux = in_compute_flux(convlik_forw[Mstep+m:Mstep+m, :, t:t],fluxt);

            // Compute metapop effects
            convlikout_forw[Mstep+m, :, t] = metapop(DO_METAPOP,DO_IN_OUT,
                forw_Rin[,m:m],forw_Rout[,m:m],convlik_forw[Mstep+m:Mstep+m, :, t:t],convforwflux,fluxproportions[m:m],fluxt)[1,:,1];

            // Draw new latent infections from observation model
            if (OBSERVATION_MODEL == POISSON) {
              Clatent[,Tcur+i] = to_vector(poisson_rng(convlikout_forw[Mstep+m, , t]));
            } else if (OBSERVATION_MODEL == NEG_BINOMIAL_2) {
              Clatent[,Tcur+i] = to_vector(neg_binomial_2_rng(convlikout_forw[Mstep+m, , t], 1 / dispersion));
            } else if (OBSERVATION_MODEL == NEG_BINOMIAL_3) {
              Clatent[,Tcur+i] = to_vector(neg_binomial_2_rng(convlikout_forw[Mstep+m, , t], convlikout_forw[Mstep+m, , t] / dispersion));
            } else if (OBSERVATION_MODEL == GAUSSIAN) {
              Clatent[,Tcur+i] = to_vector(fabs(normal_rng(convlikout_forw[Mstep+m, , t], sqrt((1.0+dispersion)*convlikout_forw[Mstep+m, , t]))));
            } else {
              reject("Invalid combination of OBSERVATION_DATA, OBSERVATION_MODEL found: ", OBSERVATION_DATA, ", ", OBSERVATION_MODEL);
            }
          }

        }
        // Simulate observations based on the type data we are observing and compute predictive probabilities
        {
          // Cases where we only observe the latent epidemic
          if (OBSERVATION_DATA == COUNTS || OBSERVATION_DATA == CLEANED_RECON){
            { // Predictions are exactly the forward simulated model
              for (t in Tcond+1:Tcur) {
                int s = t-Tcond;
                Cpred[,s] = Clatent[,t];
              }
            }
            { // Projections are exactly the forward simulated model
              for (t in Tcur+1:Tcur+Tproj) {
                int s = t-Tcur;
                Cproj[,s] = Clatent[,t];
              }
            }
            // Compute predictive likelihood of observed latent epidemic
            {
              for (m in 1:Mpred) {
                for (t in 1:Tstep) {
                  for (j in 1:N) {
                    int i = (m-1) * Tstep + t;
                    if (i > Tpred)
                      break;
                    if (OBSERVATION_DATA == COUNTS) {

                      if (OBSERVATION_MODEL == POISSON) {
                        Ppred[j,i] = exp(poisson_lpmf(Count[j,Tcur+i] |
                            convlikout_forw[m,j,t]
                        ));
                      } else if (OBSERVATION_MODEL == NEG_BINOMIAL_2) {
                        Ppred[j,i] = exp(neg_binomial_2_lpmf(Count[j,Tcur+i] |
                            convlikout_forw[m,j,t],
                            1.0 / dispersion
                        ));
                      } else if (OBSERVATION_MODEL == NEG_BINOMIAL_3) {
                        Ppred[j,i] = exp(neg_binomial_2_lpmf(Count[j,Tcur+i] |
                            convlikout_forw[m,j,t],
                            convlikout_forw[m,j,t] / dispersion
                        ));
                      } else {
                        reject("Invalid combination of OBSERVATION_DATA, OBSERVATION_MODEL found: ", OBSERVATION_DATA, ", ", OBSERVATION_MODEL);
                      }

                    } else if (OBSERVATION_DATA == CLEANED_RECON) {

                      if (OBSERVATION_MODEL == POISSON) {
                        Ppred[j,i] = exp(poisson_lpmf(Clean_recon[j,Tcur+i] |
                            convlikout_forw[m,j,t]
                        ));
                      } else if (OBSERVATION_MODEL == NEG_BINOMIAL_2) {
                        Ppred[j,i] = exp(neg_binomial_2_lpmf(Clean_recon[j,Tcur+i] |
                            convlikout_forw[m,j,t],
                            1.0 / dispersion
                        ));
                      } else if (OBSERVATION_MODEL == NEG_BINOMIAL_3) {
                        Ppred[j,i] = exp(neg_binomial_2_lpmf(Clean_recon[j,Tcur+i] |
                            convlikout_forw[m,j,t],
                            convlikout_forw[m,j,t] / dispersion
                        ));
                      } else {
                        reject("Invalid combination of OBSERVATION_DATA, OBSERVATION_MODEL found: ", OBSERVATION_DATA, ", ", OBSERVATION_MODEL);
                      }
                      
                    }
                  }
                }
              }
            }
          // Cases where we observe case reports
          } else if (OBSERVATION_DATA == CLEANED_LATENT || OBSERVATION_DATA == INFECTION_REPORTS) {
            //*** TODO Build in more observation model options, and rename OBSERVATION_MODEL to INFECTION_MODEL? ***//
            { // posterior predictive expected counts
              for (t in Tcond+1:Tcur) {
                int s = t-Tcond;
                Cpred[,s] = Clatent[,t-Tdp+1:t] * delayprofile_rev;
                // Draw sample from observation model
                Cpred[,s] = to_vector(neg_binomial_2_rng(Cpred[,s], Cpred[,s] / 1.0)); //*** TODO use better estimated dispersion ***//
              }
            }
            { // forecasting expected counts given parameters
              for (t in Tcur+1:Tcur+Tproj) {
                int s = t-Tcur;
                Cproj[,s] = Clatent[,t-Tdp+1:t] * delayprofile_rev;
                // Draw sample from observation model
                Cproj[,s] = to_vector(neg_binomial_2_rng(Cproj[,s], Cproj[,s] / 1.0)); //*** TODO use better estimated dispersion ***//
              }
            }
            // Compute predictive likelihood of observed future case observations
            {
              for (t in Tcur+1:Tcur+Tpred) {
                int i = t-Tcur;
                vector[N] cc = Clatent[,t-Tdp+1:t] * delayprofile_rev;
                for (j in 1:N) {
                  Ppred[j,i] = exp(neg_binomial_2_lpmf(Count[j,t] |
                      cc[j],
                      cc[j] / 1.0 //*** TODO use better estimated dispersion ***//
                  )); 
                }
              }
            }
          }
        }
      }
    }

    {  
      for (j in 1:N) {
        for (k in 1:(Mstep+Mproj)) {
          real s = 0.0;
          real r = 0.0;
          for (i in 1:Tstep) {
            s += convlikout_forw[k,j,i];
            r += convlik_forw[k,j,i];
          }
          convlikout_forw_reduced[k,j,1] = s + 1e-6;
          convlik_forw_reduced[k,j,1] = r + 1e-6;
        }
      }
    }

    // Estimated Rt and Rt for each and all areas
    {
      matrix[N, Mstep+Mforw] oneN = rep_matrix(1.0,N,Mstep+Mforw);
      vector[1] oneT = rep_vector(1.0,1);
      matrix[F1,N*1] convlikflux_forw_reduced[Mstep+Mforw] = in_compute_flux(convlik_forw_reduced,fluxt);
      matrix[N,1] convone_reduced[Mstep+Mforw] = metapop(DO_METAPOP,DO_IN_OUT,
          oneN,oneN,convlik_forw_reduced,convlikflux_forw_reduced,fluxproportions,fluxt);
      for (m in 1:(Mstep+Mforw)) {
        Rt_all[m] = sum(convlikout_forw_reduced[m]) / sum(convone_reduced[m]);
        Rt[m] = (convlikout_forw_reduced[m] * oneT) ./ (convone_reduced[m] * oneT);
      }
    }

  }

  

}

