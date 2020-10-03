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
  int<lower=0> Mproj;           // number of time steps
  int<lower=0> Tproj;       // number of days to forecast

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
  matrix[Mstep,Mstep] timedist;     // distance between time samples
  matrix[Mstep,Mstep] timecorcut;   // matrix specifying which time points should be correlated (to account for lockdown)

  int<lower=1,upper=5> SPATIAL_KERNEL;
  int<lower=1,upper=5> TEMPORAL_KERNEL;
  int<lower=0,upper=1> LOCAL_KERNEL;
  int<lower=0,upper=1> GLOBAL_KERNEL;
  int<lower=0,upper=1> DO_METAPOP;
  int<lower=0,upper=1> DO_IN_OUT;
  int<lower=1,upper=3> OBSERVATION_DATA;
  int<lower=1,upper=4> OBSERVATION_MODEL;
}

transformed data {
  int Tlik = Mstep*Tstep;
  int Tcur = Tcond+Tlik;    // index of day on which we are estimating Rt
  int Tpred = Tall-Tcur;    // number of days to calculate predictive probabilities for
  int F1 = F+1;
  vector[max(1,F)] ones = rep_vector(1.0,max(1,F));

  matrix[N,Tall] Creal;     // real type version of Clean
  vector[Tip] infprofile_rev; // reversed infection profile
  vector[Tdp] delayprofile_rev; // reversed infection profile
  vector[Tdp] delayprofile_cum; // reversed infection profile

  // precompute convolutions between Counts and infprofile 
  matrix[N,Tstep] convlik[Mstep];      // for use in likelihood computation
  matrix[N,1] convlik_reduced[Mstep];      // for use in likelihood computation
  matrix[N,Tpred] convpred[1];    // for use in predictive probs of future counts
  int Tforw = max(Tpred,Tproj);
  matrix[N,Tforw] convforw[1];    // for use in forecasting into future 

  matrix[F1,N*N] fluxt;      // transposed flux matrices
  matrix[F1,N*Tstep] convlikflux[Mstep];
  matrix[F1,N*1] convlikflux_reduced[Mstep];
  matrix[F1,N*Tpred] convpredflux[1];

  int Count_lik_reduced[Mstep,N];
  int Clean_recon_lik_reduced[Mstep,N];
  matrix[Mstep,N] Clean_latent_lik_reduced;

  int NONE_KERNEL = 1;
  int EXP_QUAD_KERNEL = 2;
  int MATERN12_KERNEL = 3;
  int MATERN32_KERNEL = 4;
  int MATERN52_KERNEL = 5;

  int POISSON = 1;
  int NEG_BINOMIAL_2 = 2;
  int NEG_BINOMIAL_3 = 3;
  int GAUSSIAN = 4;

  int COUNTS = 1;
  int CLEANED_LATENT = 2;
  int CLEANED_RECON = 3;

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

  for (k in 1:Mstep) {
    for (j in 1:N) {
      {
        int sum_count = 0;
        int sum_clean_recon = 0;
        real sum_clean_latent = 0;
        for (i in 1:Tstep) {
          sum_count += Count[j,Tcond+i+((k-1)*Tstep)];
          sum_clean_recon += Clean_recon[j,Tcond+i+((k-1)*Tstep)];
          sum_clean_latent += Clean_latent[j,Tcond+i+((k-1)*Tstep)];
        }
        Count_lik_reduced[k,j] = sum_count;
        Clean_recon_lik_reduced[k,j] = sum_clean_recon;
        Clean_latent_lik_reduced[k,j] = sum_clean_latent;
      }
    }
  }

  if (OBSERVATION_DATA==CLEANED_LATENT) {
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
      real s = 0.0;
      for (i in 1:Tstep) {
        int L = min(Tip,Tcond+i-1+((k-1) * Tstep)); // length of infection profile that overlaps with case counts 
        convlik[k,j,i] = dot_product(Creal[j,Tcond+i+((k-1)*Tstep)-L:Tcond+i+((k-1)*Tstep)-1], infprofile_rev[Tip-L+1:Tip]) + 1e-6;
        s += convlik[k,j,i];
      }
      convlik_reduced[k,j,1] = s + 1e-6;
    }
    
    for (i in 1:Tpred) {
      int L = min(Tip,Tcur+i-1); // length of infection profile that overlaps with case counts 
      convpred[1,j,i] = dot_product(Creal[j,Tcur-L+i:Tcur-1+i], infprofile_rev[Tip-L+1:Tip])+1e-6;
    }
    for (i in 1:Tforw) {
      int L = min(Tip,Tcur+i-1); // length of infection profile that overlaps with case counts 
      convforw[1,j,i] = dot_product(Creal[j,Tcur-L+i:Tcur], infprofile_rev[Tip-L+1:Tip-i+1])+1e-6;
    }
  }

  if (DO_METAPOP && !DO_IN_OUT) {
    convlikflux = in_compute_flux(convlik,fluxt);
    convlikflux_reduced = in_compute_flux(convlik_reduced,fluxt);
    convpredflux = in_compute_flux(convpred,fluxt);
  }
}

parameters {
  real<lower=0> gp_space_length_scale;
  real<lower=0> gp_space_sigma;

  real<lower=0> gp_time_length_scale;
    
  matrix<lower=0>[N,Mstep] local_exp;
  real<lower=0> local_scale;
  real<lower=0> global_sigma;

  vector[N*Mstep] eta_in;
  vector[N*Mstep] eta_out;

  vector[N*Mstep] epsilon_in;
  vector[N*Mstep] epsilon_out;

  real<lower=0> dispersion;
  // real<lower=0> Ravg;
  real<lower=0,upper=1> coupling_rate[Mstep];
  simplex[max(1,F)] flux_probs;
}

transformed parameters {
  matrix[N, Mstep] Rin;                 // instantaneous reproduction number
  matrix[N, Mstep] Rout;                 // instantaneous reproduction number
  matrix[N, Mstep] local_sigma;
  row_vector[F1] fluxproportions[Mstep];
  matrix[N,1] convlikout_reduced[Mstep];
  {
    matrix[N,N] K_space;
    matrix[Mstep,Mstep] K_time;
    matrix[N,N] L_space;
    matrix[Mstep,Mstep] L_time;
    real global_sigma2 = square(global_sigma);

    // GP space kernel
    K_space = kernel(SPATIAL_KERNEL,geodist, gp_space_sigma, gp_space_length_scale,
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
      for (m in 1:Mstep) {
        for (n in 1:N) {
          local_sigma[n, m] = sqrt((2 * square(local_scale)) * local_exp[n, m]);
        }
      }
    } else {
      local_sigma = rep_matrix(0.0,N,Mstep);
    }
    
//print("local_sigma");
//print(local_sigma);

    // Compute (K_time (*) K_space) * eta via efficient kronecker trick. 
    // Don't reshape back to vector for convinience.
    // Add on the location-time dependant noise as well
    Rin = exp(
        (L_space * to_matrix(eta_in, N, Mstep) * L_time') + 
        (local_sigma .* to_matrix(epsilon_in, N, Mstep))
    );
    if (DO_METAPOP && DO_IN_OUT) {
      Rout = exp(
          (L_space * to_matrix(eta_out, N, Mstep) * L_time') + 
          (local_sigma .* to_matrix(epsilon_out, N, Mstep))
      );
    } else {
      Rout = rep_matrix(1.0,N,Mstep);
    }

    // metapopulation infection rate model
    if (DO_METAPOP) {
      for (m in 1:Mstep){
        fluxproportions[m, 1] = 1.0-coupling_rate[m];
      for (f in 2:F1)
        fluxproportions[m, f] = coupling_rate[m]*flux_probs[f-1];
      }
    } else {
      for (m in 1:Mstep){
        fluxproportions[m, 1] = 1.0;
        for (f in 2:F1)
          fluxproportions[m, f] = 0.0;
      }
      
    }
    convlikout_reduced = metapop(DO_METAPOP,DO_IN_OUT,
        Rin,Rout,convlik_reduced,convlikflux_reduced,fluxproportions,fluxt);
  }
}

model {
  coupling_rate ~ normal(0.0, .1);
  flux_probs ~ dirichlet(ones);
  dispersion ~ normal(0.0,5.0);

  // GP prior density
  eta_in ~ std_normal();
  eta_out ~ std_normal();
  epsilon_in ~ std_normal();
  epsilon_out ~ std_normal();


  gp_time_length_scale ~ gig(14, 0.2, 1.0);
  gp_space_length_scale ~ gig(5, 5.0, 5.0);
  gp_space_sigma ~ normal(0.0, 0.25);
  global_sigma ~  normal(0.0, 0.25);
  local_scale ~ normal(0.0, 0.1);
  for (j in 1:Mstep){
    for (i in 1:N) {
      local_exp[i, j] ~ exponential(1.0);
      // local_sigma2[i, j] ~ exponential(0.5 / square(local_scale)); // reparameterised
    }
  }

  // compute likelihoods
  for (k in 1:Mstep) {
    for (j in 1:N) {
      if (OBSERVATION_DATA == COUNTS) {

        if (OBSERVATION_MODEL == POISSON) {
          Count_lik_reduced[k,j] ~ poisson(
              convlikout_reduced[k,j,1] 
          );
        } else if (OBSERVATION_MODEL == NEG_BINOMIAL_2) {
          Count_lik_reduced[k,j] ~ neg_binomial_2(
              convlikout_reduced[k,j,1],
              1.0 / dispersion
          );
        } else if (OBSERVATION_MODEL == NEG_BINOMIAL_3) {
          Count_lik_reduced[k,j] ~ neg_binomial_2(
              convlikout_reduced[k,j,1],
              convlikout_reduced[k,j,1] / dispersion
          );
        } else {
          reject("Invalid combination of OBSERVATION_DATA, OBSERVATION_MODEL found: ", OBSERVATION_DATA, ", ", OBSERVATION_MODEL)
        }

      } else if (OBSERVATION_DATA == CLEANED_LATENT) {

        if (OBSERVATION_MODEL == GAUSSIAN) {
          Clean_latent_lik_reduced[k,j] ~ normal(
              convlikout_reduced[k,j,1], 
              sqrt((1.0+dispersion)*convlikout_reduced[k,j,1])
          );
        } else {
          reject("Invalid combination of OBSERVATION_DATA, OBSERVATION_MODEL found: ", OBSERVATION_DATA, ", ", OBSERVATION_MODEL)
        }

      } else if (OBSERVATION_DATA == CLEANED_RECON) {

        if (OBSERVATION_MODEL == POISSON) {
          Clean_recon_lik_reduced[k,j] ~ poisson(
              convlikout_reduced[k,j,1] 
          );
        } else if (OBSERVATION_MODEL == NEG_BINOMIAL_2) {
          Clean_recon_lik_reduced[k,j] ~ neg_binomial_2(
              convlikout_reduced[k,j,1],
              1.0 / dispersion
          );
        } else if (OBSERVATION_MODEL == NEG_BINOMIAL_3) {
          Clean_recon_lik_reduced[k,j] ~ neg_binomial_2(
              convlikout_reduced[k,j,1],
              convlikout_reduced[k,j,1] / dispersion
          );
        } else {
          reject("Invalid combination of OBSERVATION_DATA, OBSERVATION_MODEL found: ", OBSERVATION_DATA, ", ", OBSERVATION_MODEL)
        }

      }
    }
  }
}

generated quantities {
  real Rt_all[Mstep];
  vector[N] Rt[Mstep];
  matrix[N,Tpred] Ppred;
  matrix[N,Mstep*Tstep] Cpred;
  matrix[N,Tproj] Cproj; 

  // Estimated Rt and Rt for each and all areas
  {
    matrix[N, Mstep] oneN = rep_matrix(1.0,N,Mstep);
    vector[1] oneT = rep_vector(1.0,1);
    matrix[N,1] convone_reduced[Mstep] = metapop(DO_METAPOP,DO_IN_OUT,
        oneN,oneN,convlik_reduced,convlikflux_reduced,fluxproportions,fluxt);
    for (m in 1:Mstep) {
      Rt_all[m] = sum(convlikout_reduced[m]) / sum(convone_reduced[m]);
      Rt[m] = (convlikout_reduced[m] * oneT) ./ (convone_reduced[m] * oneT);
    }
  }
  
  
  {
    matrix[N,Tstep] convlikout[Mstep];

    convlikout = metapop(DO_METAPOP,DO_IN_OUT,
        Rin,Rout,convlik,convlikflux,fluxproportions,fluxt);


    if (1) {
        // OBSERVATIONMODEL != CLEANED_LATENT &&
        // OBSERVATIONMODEL != CLEANED_RECON) {
      { // posterior predictive expected counts
        for (k in 1:Mstep) 
          Cpred[,(1+(k-1)*Tstep):(k*Tstep)] = convlikout[k];
      }
      { // predictive probability of future counts 
        row_vector[F1] pred_fluxproportions[1];
        matrix[N,1] pred_Rin;
        matrix[N,1] pred_Rout;
        matrix[N,Tpred] convpredout[1];

        pred_fluxproportions[1] = fluxproportions[Mstep];
        pred_Rin = block(Rin, 1, Mstep, N, 1);
        pred_Rout = block(Rout, 1, Mstep, N, 1);
        convpredout = metapop(DO_METAPOP,DO_IN_OUT,
            pred_Rin,pred_Rout,convpred,convpredflux,pred_fluxproportions,fluxt);

        for (i in 1:Tpred) {
          for (j in 1:N) {
            if (OBSERVATION_DATA == COUNTS) {

              if (OBSERVATION_MODEL == POISSON) {
                Ppred[j,i] = exp(poisson_lpmf(Count[j,Tcur+i] |
                    convpredout[1,j,i]
                ));
              } else if (OBSERVATION_MODEL == NEG_BINOMIAL_2) {
                Ppred[j,i] = exp(neg_binomial_2_lpmf(Count[j,Tcur+i] |
                    convpredout[1,j,i],
                    1.0 / dispersion
                ));
              } else if (OBSERVATION_MODEL == NEG_BINOMIAL_3) {
                Ppred[j,i] = exp(neg_binomial_2_lpmf(Count[j,Tcur+i] |
                    convpredout[1,j,i],
                    convpredout[1,j,i] / dispersion
                ));
              } else {
                reject("Invalid combination of OBSERVATION_DATA, OBSERVATION_MODEL found: ", OBSERVATION_DATA, ", ", OBSERVATION_MODEL)
              }

            } else if (OBSERVATION_DATA == CLEANED_LATENT) {

              if (OBSERVATION_MODEL == GAUSSIAN) {
                Ppred[j,i] = exp(normal_lpdf(Clean_latent[j,Tcur+i] |
                    convpredout[1,j,i],
                    sqrt((1.0+dispersion)*convpredout[1,j,i])
                ));
              } else {
                reject("Invalid combination of OBSERVATION_DATA, OBSERVATION_MODEL found: ", OBSERVATION_DATA, ", ", OBSERVATION_MODEL)
              }

            } else if (OBSERVATION_DATA == CLEANED_RECON) {

              if (OBSERVATION_MODEL == POISSON) {
                Ppred[j,i] = exp(poisson_lpmf(Clean_recon[j,Tcur+i] |
                    convpredout[1,j,i]
                ));
              } else if (OBSERVATION_MODEL == NEG_BINOMIAL_2) {
                Ppred[j,i] = exp(neg_binomial_2_lpmf(Clean_recon[j,Tcur+i] |
                    convpredout[1,j,i],
                    1.0 / dispersion
                ));
              } else if (OBSERVATION_MODEL == NEG_BINOMIAL_3) {
                Ppred[j,i] = exp(neg_binomial_2_lpmf(Clean_recon[j,Tcur+i] |
                    convpredout[1,j,i],
                    convpredout[1,j,i] / dispersion
                ));
              } else {
                reject("Invalid combination of OBSERVATION_DATA, OBSERVATION_MODEL found: ", OBSERVATION_DATA, ", ", OBSERVATION_MODEL)
              }
              
            }
          }
        }
      }
      { // forecasting expected counts given parameters
        row_vector[F1] forw_fluxproportions[1];
        matrix[F1,N*1] convforwflux[1];
        matrix[N,1] forw_Rin;
        matrix[N,1] forw_Rout;
        matrix[N,1] convforwall[1];

        forw_fluxproportions[1] = fluxproportions[Mstep];
        forw_Rin = block(Rin, 1, Mstep, N, 1);
        forw_Rout = block(Rout, 1, Mstep, N, 1);

        for (i in 1:Tproj) {
          for (j in 1:N) 
            convforwall[1,j,1] = convforw[1,j,i] + 
                dot_product(Cproj[j,1:(i-1)], infprofile_rev[(Tip-i+2):Tip]);
          if (DO_METAPOP && !DO_IN_OUT)
            convforwflux = in_compute_flux(convforwall,fluxt);
          Cproj[,i] = metapop(DO_METAPOP,DO_IN_OUT,
              forw_Rin,forw_Rout,convforwall,convforwflux,fluxproportions,fluxt)[1,:,1];
        }
      }
    } else if (0) { 
        // OBSERVATIONMODEL == CLEANED_LATENT ||
        // OBSERVATIONMODEL == CLEANED_RECON) {
      int Tidp = max(Tip,Tdp);
      matrix[N,Tcur+Tforw] Clatent;
      row_vector[F1] forw_fluxproportions[1];
      matrix[F1,N*1] convforwflux[1];
      matrix[N,1] forw_Rin;
      matrix[N,1] forw_Rout;
      matrix[N,1] convforwall[1];

      forw_fluxproportions[1] = fluxproportions[Mstep];
      forw_Rin = block(Rin, 1, Mstep, N, 1);
      forw_Rout = block(Rout, 1, Mstep, N, 1);

      Clatent[,1:Tcond] = Creal[,1:Tcond];
      for (k in 1:Mstep) 
        Clatent[,Tcond+(k-1)*Tstep+1:Tcond+k*Tstep] = convlikout[k];
      for (i in 1:Tforw) {
        for (j in 1:N) 
          convforwall[1,j,1] = convforw[1,j,i] + 
              dot_product(Clatent[j,Tcur+1:Tcur+(i-1)], infprofile_rev[(Tip-i+2):Tip]);
        if (DO_METAPOP && !DO_IN_OUT)
          convforwflux = in_compute_flux(convforwall,fluxt);
        Clatent[,Tcur+i] = metapop(DO_METAPOP,DO_IN_OUT,
            forw_Rin,forw_Rout,convforwall,convforwflux,fluxproportions,fluxt)[1,:,1];
      }
 
      { // posterior predictive expected counts
        for (t in Tcond+1:Tcur) {
          int s = t-Tcond;
          Cpred[,s] = Clatent[,t-Tdp+1:t] * delayprofile_rev;
        }
      }
      { // forecasting expected counts given parameters
        for (t in Tcur+1:Tcur+Tproj) {
          int i = t-Tcur;
          Cproj[,i] = Clatent[,t-Tdp+1:t] * delayprofile_rev;
        }
      }
      // predictive probability of future counts 
      for (t in Tcur+1:Tcur+Tpred) {
        int i = t-Tcur;
        vector[N] cc = Clatent[,t-Tdp+1:t] * delayprofile_rev;
        for (j in 1:N) {
          Ppred[j,i] = exp(neg_binomial_2_lpmf(Count[j,t] |
              cc[j],
              cc[j] / 2.0 //*** TODO use better estimated dispersion ***//
          )); 
        }
      }
    }

  } 

}

