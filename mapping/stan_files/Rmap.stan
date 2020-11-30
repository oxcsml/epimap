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
  int<lower=1> Tall;        // number of all days in case count time series
  int<lower=1> Tcur;       // number of days we will condition on
  int<lower=0> Tstep;       // number of days to step for each time step of Rt prediction
  int<lower=0> Tpred;       // number of days to compute predictive likelihoods for
  int<lower=1> Mstep;           // number of time steps
  int<lower=0> Mignore;           // number of time steps to ignore
  int<lower=0> Mproj;           // number of time steps

  int Ct[N,Tall];           // case counts
  int Clean_recon[N,Tcur];     // case counts
  matrix[N,Tcur] Clean_latent; // cleaned case counts

  // vector[2] geoloc[N];   // geo locations of regions
  int<lower=1> Tip;         // length of infection profile
  vector[Tip] infprofile;   // infection profile aka serial interval distribution
  int<lower=1> Tdp;         // length of infection profile
  vector[Tdp] delayprofile; // infection profile aka serial interval distribution
  int F;
  matrix[N,N] flux[F];      // fluxes for radiation metapopulation model

  int<lower=1> N_region;               // number of regions
  matrix[N,N_region] sparse_region; 

  matrix[N,N] geodist;      // distance between locations
  matrix[Mstep+Mproj,Mstep+Mproj] timedist;     // distance between time samples
  matrix[Mstep+Mproj,Mstep+Mproj] timecorcut;   // matrix specifying which time points should be correlated (to account for lockdown)

  real<lower=0> gp_space_scale; // 50km
  real<lower=0> gp_time_scale; // 3 weeks
  real<lower=0> gp_space_decay_scale;
  real<lower=0> gp_time_decay_scale;
  real fixed_gp_space_length_scale; // If positive, use this, otherwise use prior
  real fixed_gp_time_length_scale; // If positive, use this, otherwise use prior

  real coupling_mu_loc;
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
  int<lower=0,upper=1> CONSTANT_FORWARD_RT;
  int<lower=0,upper=1> FULL_CASES_DISTRIBUTION;
}

transformed data {
  int Tlik = Mstep*Tstep;     // Number of days to use for the likelihood computation
  int Tproj = Mproj*Tstep;    // The number of days to project forward
  int Tcond = Tcur-Tlik;    // index of day on which we are estimating Rt
  int Mpred = (Tpred + Tstep - 1) / Tstep; // Enough Tsteps to cover whole Tpred period. Funky math to deal with Stans deficiencies
  int Tforw = max(Tpred,Tproj);
  int Mforw = max(Mpred,Mproj);
  int F1 = F+1;
  vector[max(1,F)] ones = rep_vector(1.0,max(1,F));

  matrix[N,Tcur] Xt;     // real type version of Clean
  vector[Tip] infprofile_rev; // reversed infection profile
  vector[Tdp] delayprofile_rev; // reversed infection profile
  vector[Tdp] delayprofile_cum; // reversed infection profile

  // precompute convolutions between Counts and infprofile 
  matrix[N,Tstep] Zt[Mstep];      // for use in likelihood computation
  // YW COMMENTED OUT
  //matrix[N,Tforw] convforw[1];    // for use in forecasting into future 

  matrix[F1,N*N] fluxt;      // transposed flux matrices
  matrix[F1,N*Tstep] FZt[Mstep];

  matrix[N,N] fixed_L_space;
  matrix[Mstep+Mforw,Mstep+Mforw] fixed_L_time;

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
    Xt = Clean_latent;
  } else if (OBSERVATION_DATA==CLEANED_RECON) {
    for (j in 1:N) 
      for (i in 1:Tcur)
        Xt[j,i] = Clean_recon[j,i];
  } else {
    for (j in 1:N) 
      for (i in 1:Tcur)
        Xt[j,i] = Ct[j,i];
  }

  // precompute convolutions between counts and infprofile
  // compute for each time offset - NOTE: not the fastest ordering of these access options (see https://mc-stan.org/docs/2_23/stan-users-guide/indexing-efficiency-section.html), but assuming okay as this is only done once
  for (j in 1:N) {
    for (k in 1:Mstep) {
      for (i in 1:Tstep) {
        int L = min(Tip,Tcond+i-1+((k-1) * Tstep)); // length of infection profile that overlaps with case counts 
        Zt[k,j,i] = dot_product(Xt[j,Tcond+i+((k-1)*Tstep)-L:Tcond+i+((k-1)*Tstep)-1], infprofile_rev[Tip-L+1:Tip]) + 1e-6;
      }
    }
    
    //for (i in 1:Tpred) {
    //  int L = min(Tip,Tcur+i-1); // length of infection profile that overlaps with case counts 
    //  convpred[1,j,i] = dot_product(Xt[j,Tcur-L+i:Tcur-1+i], infprofile_rev[Tip-L+1:Tip])+1e-6;
    //}
    // YW COMMENTED OUT
    //for (i in 1:Tforw) {
    //  int L = min(Tip,Tcur+i-1); // length of infection profile that overlaps with case counts 
    //  convforw[1,j,i] = dot_product(Xt[j,Tcur-L+i:Tcur], infprofile_rev[Tip-L+1:Tip-i+1])+1e-6;
    //}
  }

  if (DO_METAPOP && !DO_IN_OUT) {
    FZt = in_compute_flux(Zt,fluxt);
    //convpredflux = in_compute_flux(convpred,fluxt);
  }

  if (fixed_gp_space_length_scale > 0.0) {
    matrix[N,N] fixed_K_space;
    fixed_K_space = kernel(SPATIAL_KERNEL,geodist, 1.0, fixed_gp_space_length_scale,
        NONE_KERNEL,EXP_QUAD_KERNEL,MATERN12_KERNEL,MATERN32_KERNEL,MATERN52_KERNEL);
    fixed_L_space = cholesky_decompose(fixed_K_space);
  }

  if (fixed_gp_time_length_scale > 0.0) {
    matrix[Mstep+Mforw,Mstep+Mforw] fixed_K_time;
    fixed_K_time  = kernel(TEMPORAL_KERNEL,timedist, 1.0, fixed_gp_time_length_scale,
        NONE_KERNEL,EXP_QUAD_KERNEL,MATERN12_KERNEL,MATERN32_KERNEL,MATERN52_KERNEL);
    fixed_K_time  = fixed_K_time .* timecorcut;  // Zero out uncorrelated time entries
    fixed_L_time = cholesky_decompose(fixed_K_time)';
  }

}

parameters {
  real<lower=0.0> gp_sigma;
  real<lower=0.0,upper=0.632> gp_space_decay;
  real<lower=0.0,upper=0.632> gp_time_decay;
    
  matrix<lower=0.0>[LOCAL_KERNEL ? N : 1,LOCAL_KERNEL ? Mstep+Mforw : 1] local_exp;
  real<lower=0.0> local_scale;
  real<lower=0.0> global_sigma;

  row_vector[(Mstep+Mforw)] global_eta_in;
  row_vector[(DO_METAPOP && DO_IN_OUT) ? (Mstep+Mforw) : 1] global_eta_out;

  vector[N*(Mstep+Mforw)] gp_eta_in;
  vector[(DO_METAPOP && DO_IN_OUT) ? N*(Mstep+Mforw) : 1] gp_eta_out;

  vector[N*(Mstep+Mforw)] local_eta_in;
  vector[(DO_METAPOP && DO_IN_OUT) ? N*(Mstep+Mforw) : 1] local_eta_out;

  real<lower=0.0> infection_dispersion;
  real<lower=0.0> case_dispersion;

  real coupling_mu_eta; // prior mean for sigmoid^-1(coupling_rate) 
  real<lower=0> coupling_sigma; // prior scale for sigmoid^-1(coupling_rate) 
  real<lower=0,upper=1> coupling_alpha1; // 1-autocorrelation for sigmoid^-1(coupling_rate) 
  vector[DO_METAPOP ? (Mstep+Mforw) : 1] coupling_eta; // noise for sigmoid^-1(coupling_rate) 

  simplex[max(1,F)] flux_probs;
}

transformed parameters {
  matrix[N, Mstep+Mforw] Rin;                 // instantaneous reproduction number
  matrix[N, Mstep+Mforw] Rout;                 // instantaneous reproduction number
  row_vector[F1] fluxproportions[Mstep+Mforw];
  matrix[N,Tstep] EXt[Mstep];
  real gp_space_length_scale;
  real gp_time_length_scale;
  real coupling_alpha = 1.0-coupling_alpha1; // autocorrelation for sigmoid^-1(coupling_rate) 
  real coupling_rate[Mstep+Mforw];

  if (DO_METAPOP) {
    // AR1 process for sigmoid^-1(coupling_rate)
    vector[Mstep+Mforw] Yt = rep_vector(0.0,Mstep+Mforw);
    real delta = sqrt(1.0-pow(coupling_alpha, 2));
    Yt[1] += coupling_eta[1];
    for (i in 2:(Mstep+Mforw))
      Yt[i] += coupling_alpha * Yt[i-1] + delta * coupling_eta[i];
    coupling_rate = to_array_1d(inv_logit(
      coupling_mu_loc + 
      coupling_mu_scale * coupling_mu_eta + 
      coupling_sigma * Yt
    ));
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

  { // Prior on Rt
  matrix[N,N] L_space;
  matrix[Mstep+Mforw,Mstep+Mforw] L_time;
  if (fixed_gp_space_length_scale <= 0.0) {
    matrix[N,N] K_space;

    gp_space_length_scale = - gp_space_scale / log(1.0-gp_space_decay);
    
    K_space = kernel(SPATIAL_KERNEL,geodist, 1.0, gp_space_length_scale,
        NONE_KERNEL,EXP_QUAD_KERNEL,MATERN12_KERNEL,MATERN32_KERNEL,MATERN52_KERNEL);
    L_space = cholesky_decompose(K_space);
  } else {
    gp_space_length_scale = fixed_gp_space_length_scale;
    L_space = fixed_L_space;
  }

  if (fixed_gp_time_length_scale <= 0.0) {
    matrix[Mstep+Mforw,Mstep+Mforw] K_time;

    gp_time_length_scale = - gp_time_scale / log(1.0-gp_time_decay);

    K_time  = kernel(TEMPORAL_KERNEL,timedist, 1.0, gp_time_length_scale,
        NONE_KERNEL,EXP_QUAD_KERNEL,MATERN12_KERNEL,MATERN32_KERNEL,MATERN52_KERNEL);
    K_time  = K_time .* timecorcut;  // Zero out uncorrelated time entries

    L_time = cholesky_decompose(K_time)';
  } else {
    gp_time_length_scale  = fixed_gp_time_length_scale;
    L_time  = fixed_L_time;
  }

  {
    matrix[N, Mstep+Mforw] global_effect_in;
    matrix[N, Mstep+Mforw] local_effect_in;
    // Noise kernel (reshaped from (N*Mstep, N*Mstep) to (N,Mstep) since diagonal)
    if (GLOBAL_KERNEL) {
      global_effect_in = rep_matrix(global_sigma * global_eta_in * L_time, N);
    } else {
      global_effect_in = rep_matrix(0.0,N,Mstep+Mforw);
    }
    if (LOCAL_KERNEL) {
      matrix[N, Mstep+Mforw] local_sigma = sqrt((2.0 * square(local_scale)) * local_exp);
      local_effect_in = local_sigma .* to_matrix(local_eta_in, N, Mstep+Mforw);
    } else {
      local_effect_in = rep_matrix(0.0,N,Mstep+Mforw);
    }

    // Compute (K_time (*) K_space) * eta via efficient kronecker trick. 
    // Don't reshape back to vector for convinience.
    // Add on the location-time dependant noise as well
    if (0) { // AR1 process for GP temporal structure
      real alpha = exp(-Tstep/gp_time_length_scale);
      real delta = sqrt(1.0-alpha*alpha);
      matrix[N,Mstep+Mforw] gp_eta = to_matrix(gp_eta_in,N,Mstep+Mforw);
      matrix[N,Mstep+Mforw] gp_Yt;
      row_vector[Mstep+Mforw] global_Yt;

      gp_Yt[,1] = gp_eta[,1];
      for (i in 2:(Mstep+Mforw))
        gp_Yt[,i] = alpha * gp_Yt[,i-1] + delta * gp_eta[,i];

      global_Yt[1] = global_eta_in[1];
      for (i in 2:(Mstep+Mforw))
        global_Yt[i] = alpha * global_Yt[i-1] + delta * global_eta_in[i];
      Rin = exp(
        (gp_sigma * L_space * gp_Yt) 
        + rep_matrix(global_sigma * global_Yt, N)
        + local_effect_in
      );
      if (0) {
        print("Rin ", 
          max(Rin[,Mstep-1]), ", ", 
          max(Rin[,Mstep]), ", ", 
          max(Rin[,Mstep+Mforw])
        );
        print("gp_Yt ", 
          max(gp_Yt[,Mstep-1]), ", ", 
          max(gp_Yt[,Mstep]), ", ", 
          max(gp_Yt[,Mstep+Mforw])
        );
        print("global_Yt ", 
          (global_Yt[Mstep-1]), ", ", 
          (global_Yt[Mstep]), ", ", 
          (global_Yt[Mstep+Mforw])
        );
        print("local_Yt ", 
          max(local_effect_in[,Mstep-1]), ", ", 
          max(local_effect_in[,Mstep]), ", ", 
          max(local_effect_in[,Mstep+Mforw])
        );
      }
    } else { // Kronecker decomposition
      Rin = exp(
        (gp_sigma * L_space * to_matrix(gp_eta_in, N, Mstep+Mforw) * L_time) 
        + global_effect_in 
        + local_effect_in
      );
    }

    if (DO_METAPOP && DO_IN_OUT) {
      matrix[N, Mstep+Mforw] global_effect_out;
      matrix[N, Mstep+Mforw] local_effect_out;
      // Noise kernel (reshaped from (N*Mstep, N*Mstep) to (N,Mstep) since diagonal)
      if (GLOBAL_KERNEL) {
        global_effect_out = rep_matrix(global_sigma * global_eta_out * L_time, N);
      } else {
        global_effect_out = rep_matrix(0.0,N,Mstep+Mforw);
      }
      if (LOCAL_KERNEL) {
        matrix[N, Mstep+Mforw] local_sigma = sqrt((2.0 * square(local_scale)) * local_exp);
        local_effect_out = local_sigma .* to_matrix(local_eta_out, N, Mstep+Mforw);
      } else {
        local_effect_out = rep_matrix(0.0,N,Mstep+Mforw);
      }
      Rout = exp(
        (gp_sigma * L_space * to_matrix(gp_eta_out, N, Mstep+Mforw) * L_time) 
        + global_effect_out 
        + local_effect_out
      );
    } else {
      Rout = rep_matrix(1.0,N,Mstep+Mforw);
    }

  }
  }

  if (CONSTANT_FORWARD_RT) {
    for (m in (Mstep+1):(Mstep+Mforw)) {
      Rin[:,m] = Rin[,Mstep];
      Rout[:,m] = Rout[,Mstep];
    }
  }

    // metapopulation infection rate model
  EXt = metapop(DO_METAPOP,DO_IN_OUT,
      block(Rin, 1, 1, N, Mstep),block(Rout, 1, 1, N, Mstep),Zt,FZt,fluxproportions[1:Mstep],fluxt);
    // print("Testing infered Rts: ")
    // for (m in 1:Mstep){
    //   for (i in 1:N) {
    //     if (Rin[i, m] > 10){
    //       print("High Rt at ", i, " ", m, " value: ", Rin[i, m]);
    //       print("Rt for region: ", Rin[i]);
    //       print("Convlikout: ", EXt[m, i]);
    //       print("Global effect: ", exp(global_effect_in[i]))
    //       print("Local effect: ", exp(local_effect_in[i]))
    //     }
    //   }
    // }
    // print("Testing forward Rts: ")
    // for (m in (Mstep+1):(Mstep+Mforw)){
    //   for (i in 1:N) {
    //     if (Rin[i, m] > 10){
    //       print("High Rt at ", i, " ", m, " value: ", Rin[i, m]);
    //       print("Rt for region: ", Rin[i]);
    //       // print("Convlikout: ", EXt[m, i]);
    //       print("Global effect: ", exp(global_effect_in[i]))
    //       print("Local effect: ", exp(local_effect_in[i]))
    //     }
    //   }
    // }
}

model {
  flux_probs ~ dirichlet(ones);
  infection_dispersion ~ normal(0.0,5.0);
  case_dispersion ~ normal(0.0,5.0);

  // GP prior density
  gp_eta_in ~ std_normal();
  gp_eta_out ~ std_normal();
  local_eta_in ~ std_normal();
  local_eta_out ~ std_normal();
  global_eta_in ~ std_normal();
  global_eta_out ~ std_normal();


  gp_space_decay ~ normal(0.0,gp_space_decay_scale);
  gp_time_decay ~ normal(0.0,gp_time_decay_scale);
  gp_sigma ~ normal(0.0, 0.5);
  global_sigma ~  normal(0.0, 0.5);
  local_scale ~ normal(0.0, 0.2);
  if (LOCAL_KERNEL) {
    for (j in 1:(Mstep+Mforw))
      for (i in 1:N) 
        local_exp[i, j] ~ exponential(1.0);
  } else {
    local_exp[1,1] ~ exponential(1.0);
  }
  

  
  coupling_mu_eta ~ std_normal();
  coupling_sigma ~ normal(0.0, coupling_sigma_scale);
  coupling_alpha1 ~ normal(0.0, coupling_alpha_scale);
  coupling_eta ~ std_normal();

  // compute likelihoods
  // TODO INFECTION_REPORTS and CLEANED_LATENT code below can be merged to mostly not be duplicated
  if (OBSERVATION_DATA == INFECTION_REPORTS) {
    matrix[N,Tcur] Elatent;
    Elatent[,1:Tcond] = Xt[,1:Tcond];
    for (k in 1:Mstep) {
      Elatent[,(Tcond+1+(k-1)*Tstep):(Tcond+k*Tstep)] = EXt[k];
      if (0 && k <= Mstep-Mignore) { // TODO don't use this likelihood
        for (i in 1:Tstep) {
          int t = Tcond+i+(k-1)*Tstep;
          for (j in 1:N) {
            real Xlatent = Xt[j,t];
            real EXlatent = Elatent[j,t];
            target += log(1.0+exp(-2.0*Xlatent/(1.0+infection_dispersion))) +
                normal_lpdf(Xlatent|EXlatent,sqrt((1.0+infection_dispersion)*EXlatent));
          }
        }
      }
    }
    for (s in 1:Tlik) {
      int t = Tcond+s;
      vector[N] Ereport = Elatent[,t-Tdp+1:t] * delayprofile_rev;
      for (j in 1:N) {
        Ct[j,t] ~ neg_binomial_2(
              Ereport[j],
              Ereport[j] / case_dispersion
        ); 
      }
    }
  } else if (OBSERVATION_DATA == CLEANED_LATENT) {
    if (OBSERVATION_MODEL == GAUSSIAN) {
      for (k in 1:Mstep-Mignore) {
        for (i in 1:Tstep) {
          int t = Tcond+i+(k-1)*Tstep;
          for (j in 1:N) {
            real Xlatent = Xt[j,t];
            real EXlatent = EXt[k,j,i];
            //real loglik1;
            //real loglik2;
            //if (0) {
            //  vector[2] loglik;
            //  loglik[1] = normal_lpdf( Xlatent | Elatent, sqrt((1.0+infection_dispersion)*Elatent) );
            //  loglik[2] = normal_lpdf( -Xlatent | Elatent, sqrt((1.0+infection_dispersion)*Elatent) );
            //  loglik1 = log_sum_exp(loglik);
            //}
            //if (1) { // below code may be faster?
            //  loglik2 = log(1.0+exp(-2.0*Xlatent/(1.0+infection_dispersion))) +
            //      normal_lpdf(Xlatent|EXlatent,sqrt((1.0+infection_dispersion)*EXlatent));
            //}
            //if (fabs(loglik1-loglik2)>1e-10)
            //  print(fabs(loglik1-loglik2));
            target += log(1.0+exp(-2.0*Xlatent/(1.0+infection_dispersion))) +
                normal_lpdf(Xlatent|EXlatent,sqrt((1.0+infection_dispersion)*EXlatent));

          }
        }
      }
      for (s in 1:Tlik) {
        int t = Tcond+s;
        vector[N] Ereport = Xt[,t-Tdp+1:t] * delayprofile_rev;
        for (j in 1:N) {
          Ct[j,t] ~ neg_binomial_2(
                Ereport[j],
                Ereport[j] / case_dispersion
          );
        }
      }

          // Clean_latent_lik_reduced[k,j] ~ normal(
          //     EXt_reduced[k,j,1], 
          //     sqrt((1.0+infection_dispersion)*EXt_reduced[k,j,1])
          // );
    } else {
      reject(
        "Invalid combination of OBSERVATION_DATA, OBSERVATION_MODEL found: ", 
        OBSERVATION_DATA, ", ", OBSERVATION_MODEL
      );
    }
  } else if (OBSERVATION_DATA == COUNTS) {
    if (OBSERVATION_MODEL == POISSON) {
      for (k in 1:Mstep) {
        // Compute infinatly divisible weekly likelihood
        vector[N] EXt_weekly = rep_vector(0.0, N);
        int weekly_cases[N] = rep_array(0, N);

        for (s in 1:Tstep) {
          int t = Tcond + Tstep*(k-1)+s;
          EXt_weekly += Xt[,t-Tdp+1:t] * delayprofile_rev;
          for (j in 1:N) {
            weekly_cases[j] += Ct[j,t];
          }
        }

        for (j in 1:N) {
          weekly_cases[j] ~ poisson(
              EXt_weekly[j]
          );
        }
      }
    } else {
      reject("Currently invalid - needs refactor to daily likelihoods");
    }
    // } else if (OBSERVATION_MODEL == NEG_BINOMIAL_2) {
    //   for (k in 1:Mstep) 
    //     for (j in 1:N) 
    //       Ct[k,j] ~ neg_binomial_2(
    //           EXt_reduced[k,j,1],
    //           1.0 / case_dispersion
    //       );
    // } else if (OBSERVATION_MODEL == NEG_BINOMIAL_3) {
    //   for (k in 1:Mstep) 
    //     for (j in 1:N) 
    //       Ct[k,j] ~ neg_binomial_2(
    //           EXt_reduced[k,j,1],
    //           EXt_reduced[k,j,1] / case_dispersion
    //       );
    // } else {
    //   reject(
    //     "Invalid combination of OBSERVATION_DATA, OBSERVATION_MODEL found: ", 
    //     OBSERVATION_DATA, ", ", OBSERVATION_MODEL
    //   );
    // }

  } else if (OBSERVATION_DATA == CLEANED_RECON) {
    reject("Currently invalid - needs refactor to daily likelihoods");
    // if (OBSERVATION_MODEL == POISSON) {
    //   for (k in 1:Mstep) 
    //     for (j in 1:N) 
    //       Clean_recon_lik_reduced[k,j] ~ poisson(
    //           EXt_reduced[k,j,1] 
    //       );
    // } else if (OBSERVATION_MODEL == NEG_BINOMIAL_2) {
    //   for (k in 1:Mstep) 
    //     for (j in 1:N) 
    //       Clean_recon_lik_reduced[k,j] ~ neg_binomial_2(
    //           EXt_reduced[k,j,1],
    //           1.0 / case_dispersion
    //       );
    // } else if (OBSERVATION_MODEL == NEG_BINOMIAL_3) {
    //   for (k in 1:Mstep) 
    //     for (j in 1:N) 
    //       Clean_recon_lik_reduced[k,j] ~ neg_binomial_2(
    //           EXt_reduced[k,j,1],
    //           EXt_reduced[k,j,1] / case_dispersion
    //       );
    // } else {
    //   reject(
    //     "Invalid combination of OBSERVATION_DATA, OBSERVATION_MODEL found: ", 
    //     OBSERVATION_DATA, ", ", OBSERVATION_MODEL
    //   );
    // }
  }
}

generated quantities {
  real Rt_all[Mstep+Mproj];
  vector[N] Rt[Mstep+Mproj];
  vector[N_region] Rt_region[Mstep+Mproj];
  matrix[N,Tpred] Ppred;
  matrix[N,Mstep*Tstep] Cpred;
  matrix[N,Mproj*Tstep] Cproj;
  matrix[N_region,Mstep*Tstep] Cpred_region;
  matrix[N_region,Mproj*Tstep] Cproj_region;

  {
    matrix[N,Tcur+Tforw] Clatent;     // Latent epidemic values
    matrix[N,Tstep] Zt_forw[Mstep+Mforw]; // Infection pressure
    matrix[N,1] Zt_forw_reduced[Mstep+Mforw]; // for use in Rt computation
    matrix[N,Tstep] EXt_forw[Mstep+Mforw]; // Infection pressure
    matrix[N,1] EXt_forw_reduced[Mstep+Mforw]; // for use in Rt computation

    // Fill Clatent with observations up until the current time
    Clatent[,1:Tcur] = Xt[,1:Tcur];
    for (k in 1:Mstep) {
      EXt_forw[k] = EXt[k];
      Zt_forw[k] = Zt[k];
    }

    // Forward simulate model and compute predictive probabilities
    // Rollout stochasitc prediction of epidemic
    if (Mforw > 0) {
      // Rollout stochasitc prediction of epidemic
      {
        matrix[F1,N*1] FZt_forw[1];
        matrix[N,Mforw] forw_Rin;
        matrix[N,Mforw] forw_Rout;        

        // Pick out the posterior samples of the future parameters conditioned on past observations
        forw_Rin = block(Rin, 1, Mstep+1, N, Mforw);
        forw_Rout = block(Rout, 1, Mstep+1, N, Mforw);

        // For each timestep simulate the model forwards
        // print("Testing rollout Rts: ")
        for (m in 1:Mforw) {
          for (t in 1:Tstep) {
            int i = (m-1) * Tstep + t;
            // Compute the infection pressure by adding on pressure from simulated counts
            //print("forw_Rin ",max(forw_Rin[,m]));
            //print("fluxproportions ",fluxproportions[m]);
            for (j in 1:N) 
              Zt_forw[Mstep+m, j, t] = // convforw[1,j,i] + // YW EDITED
                  dot_product(Clatent[j,Tcur+i-Tip:Tcur+i-1], infprofile_rev);
            //print("Zt_forw ",max(Zt_forw[Mstep+m,,t]));

            // Compute flux
            if (DO_METAPOP && !DO_IN_OUT)
              FZt_forw = in_compute_flux(Zt_forw[Mstep+m:Mstep+m, :, t:t],fluxt);
            //print("convforwflux ",max(convforwflux[1,1,]));

            // Compute metapop effects
            EXt_forw[Mstep+m, :, t] = metapop(
              DO_METAPOP,DO_IN_OUT,
              forw_Rin[,m:m],
              forw_Rout[,m:m],
              Zt_forw[Mstep+m:Mstep+m, :, t:t],
              FZt_forw,
              fluxproportions[m:m], 
              fluxt
            )[1,:,1];
            //print("EXt_forw ",max(EXt_forw[Mstep+m,,t]));

            // for (n in 1:N) {
            //   if (forw_Rin[n, m] > 10){
            //     print("High forward Rt at ", n, " ", m, " value: ", forw_Rin[n, m])
            //     print("Rt for region: ", Rin[n])
            //     print("Convlikout: ", EXt_forw[Mstep+m, n])
            //   }
            // }
            // Draw new latent infections from observation model
            if (OBSERVATION_DATA == INFECTION_REPORTS || OBSERVATION_DATA == CLEANED_LATENT) {
              // TODO
              real dispersion = OBSERVATION_DATA == INFECTION_REPORTS ? 2.35 : infection_dispersion;
              Clatent[,Tcur+i] = to_vector(fabs(
                EXt_forw[Mstep+m, , t] +
                to_vector(normal_rng(rep_vector(0.0,N),rep_vector(1.0,N))) .*
                sqrt((1.0+dispersion) * EXt_forw[Mstep+m, , t])
              ));
              //print("Clatent ",Clatent[1,Tcur+i]);
            // else OBSERVATION_DATA is COUNTS or CLEANED_RECON
            } else if (OBSERVATION_MODEL == POISSON) {
              Clatent[,Tcur+i] = to_vector(poisson_rng(EXt_forw[Mstep+m, , t]));
            } else if (OBSERVATION_MODEL == NEG_BINOMIAL_2) {
              Clatent[,Tcur+i] = to_vector(neg_binomial_2_rng(EXt_forw[Mstep+m, , t], 1.0 / case_dispersion));
            } else if (OBSERVATION_MODEL == NEG_BINOMIAL_3) {
              Clatent[,Tcur+i] = to_vector(neg_binomial_2_rng(EXt_forw[Mstep+m, , t], EXt_forw[Mstep+m, , t] / case_dispersion));
            } else {
              reject("Invalid combination of OBSERVATION_DATA, OBSERVATION_MODEL found: ", OBSERVATION_DATA, ", ", OBSERVATION_MODEL);
            }
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
              for (n in 1:N_region)
                Cpred_region[n,s] = sum(Cpred[,s] .* sparse_region[,n]);
            }
          }
          { // Projections are exactly the forward simulated model
            for (t in Tcur+1:Tcur+Tproj) {
              int s = t-Tcur;
              Cproj[,s] = Clatent[,t];
              for (n in 1:N_region)
                Cproj_region[n,s] = sum(Cproj[,s] .* sparse_region[,n]);
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
                        Ppred[j,i] = exp(poisson_lpmf(Ct[j,Tcur+i] |
                            EXt_forw[m,j,t]
                        ));
                      } else if (OBSERVATION_MODEL == NEG_BINOMIAL_2) {
                        Ppred[j,i] = exp(neg_binomial_2_lpmf(Ct[j,Tcur+i] |
                            EXt_forw[m,j,t],
                            1.0 / case_dispersion
                        ));
                      } else if (OBSERVATION_MODEL == NEG_BINOMIAL_3) {
                        Ppred[j,i] = exp(neg_binomial_2_lpmf(Ct[j,Tcur+i] |
                            EXt_forw[m,j,t],
                            EXt_forw[m,j,t] / case_dispersion
                        ));
                      } else {
                        reject("Invalid combination of OBSERVATION_DATA, OBSERVATION_MODEL found: ", OBSERVATION_DATA, ", ", OBSERVATION_MODEL);
                      }

                    } else if (OBSERVATION_DATA == CLEANED_RECON) {

                      if (OBSERVATION_MODEL == POISSON) {
                        Ppred[j,i] = exp(poisson_lpmf(Clean_recon[j,Tcur+i] |
                            EXt_forw[m,j,t]
                        ));
                      } else if (OBSERVATION_MODEL == NEG_BINOMIAL_2) {
                        Ppred[j,i] = exp(neg_binomial_2_lpmf(Clean_recon[j,Tcur+i] |
                            EXt_forw[m,j,t],
                            1.0 / case_dispersion
                        ));
                      } else if (OBSERVATION_MODEL == NEG_BINOMIAL_3) {
                        Ppred[j,i] = exp(neg_binomial_2_lpmf(Clean_recon[j,Tcur+i] |
                            EXt_forw[m,j,t],
                            EXt_forw[m,j,t] / case_dispersion
                        ));
                      } else {
                        reject("Invalid combination of OBSERVATION_DATA, OBSERVATION_MODEL found: ", OBSERVATION_DATA, ", ", OBSERVATION_MODEL);
                      }
                      
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
              if (FULL_CASES_DISTRIBUTION) {
                Cpred[,s] = to_vector(neg_binomial_2_rng(Cpred[,s], Cpred[,s] / case_dispersion)); //*** TODO use better estimated dispersion ***//
              }
              for (n in 1:N_region)
                Cpred_region[n,s] = sum(Cpred[,s] .* sparse_region[,n]);
            }
          }
          { // forecasting expected counts given parameters
            for (t in Tcur+1:Tcur+Tproj) {
              int s = t-Tcur;
              Cproj[,s] = Clatent[,t-Tdp+1:t] * delayprofile_rev;
              // Draw sample from observation model
              if (FULL_CASES_DISTRIBUTION) {
                Cproj[,s] = to_vector(neg_binomial_2_rng(Cproj[,s], Cproj[,s] / case_dispersion)); //*** TODO use better estimated dispersion ***//
              }
              for (n in 1:N_region)
                Cproj_region[n,s] = sum(Cproj[,s] .* sparse_region[,n]);
            }
          }
          // Compute predictive likelihood of observed future case observations
          {
            for (t in Tcur+1:Tcur+Tpred) {
              int i = t-Tcur;
              vector[N] cc = Clatent[,t-Tdp+1:t] * delayprofile_rev;
              for (j in 1:N) {
                Ppred[j,i] = exp(neg_binomial_2_lpmf(Ct[j,t] |
                    cc[j],
                    cc[j] / case_dispersion //*** TODO use better estimated dispersion ***//
                )); 
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
            s += EXt_forw[k,j,i];
            r += Zt_forw[k,j,i];
          }
          EXt_forw_reduced[k,j,1] = s + 1e-6;
          Zt_forw_reduced[k,j,1] = r + 1e-6;
        }
      }
    }

    // Estimated Rt and Rt for each and all areas
    {
      matrix[N, Mstep+Mproj] oneN = rep_matrix(1.0,N,Mstep+Mforw);
      vector[1] oneT = rep_vector(1.0,1);
      matrix[F1,N*1] FZt_forw_reduced[Mstep+Mforw] = in_compute_flux(Zt_forw_reduced,fluxt);
      matrix[N,1] EXt_one_reduced[Mstep+Mforw] = metapop(DO_METAPOP,DO_IN_OUT,
          oneN,oneN,Zt_forw_reduced,FZt_forw_reduced,fluxproportions,fluxt);
      for (m in 1:(Mstep+Mforw)) {
        Rt_all[m] = sum(EXt_forw_reduced[m]) / sum(EXt_one_reduced[m]);
        Rt[m] = (EXt_forw_reduced[m] * oneT) ./ (EXt_one_reduced[m] * oneT);
        for (n in 1:N_region) {
          matrix[N,1] region_slice = block(sparse_region, 1, n, N, 1);
          Rt_region[m,n] = sum(EXt_forw_reduced[m] .* region_slice) / sum(EXt_forw_reduced[m] .* region_slice);
        }
      }
    }

  }

  { // print stats
    if (uniform_rng(0.0,1.0)<0.1) {
      print(
        "space ", gp_space_length_scale,
        "; time ", gp_time_length_scale,
        "; sigmas gp ", gp_sigma," g ", global_sigma," l ", local_scale,
        "; dispersion X ",infection_dispersion," C ",case_dispersion,
        "; last Rt ", Rt_all[Mstep-1], ", ", Rt_all[Mstep], ", ", Rt_all[Mstep+Mforw],
        "; last Rin ", Rin[1,Mstep-1], ", ", Rin[1,Mstep], ", ", Rin[1,Mstep+Mforw]
      );
    }
  }
}

