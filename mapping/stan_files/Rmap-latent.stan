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
    return diag_matrix(rep_vector(1e-6, cols(dist))); // very small diagonal to allow cholesky factorisation
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

  vector[Tip] infprofile_rev; // reversed infection profile
  vector[Tdp] delayprofile_rev; // reversed infection profile
  int mean_serial_interval_int = 0;
  real mean_serial_interval_real;
  matrix[N,Tcond] Xcond;

  matrix[F1,N*N] fluxt;      // transposed flux matrices

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
    real mean_serial_interval = 0.0;
    for (i in 1:Tip) {
      infprofile_rev[i] = infprofile[Tip-i+1];
      mean_serial_interval += i*infprofile[i];
    }
    // find the integer component and remainder of mean_serial_interval
    for (i in 1:Tip) {
      if (mean_serial_interval > 1.0-1e-10) {
        mean_serial_interval_int += 1;
        mean_serial_interval -= 1.0;
      }
    }
    mean_serial_interval_real = fmax(0.0,mean_serial_interval);
    for (i in 1:Tdp) {
      delayprofile_rev[Tdp-i+1] = delayprofile[i];
    }
  }

  { // mean shift the counts for rough estimate of incidence rates in past
    for (i in 1:N) {
      for (t in 1:Tcond) {
        Xcond[i,t] = (
          (1.0-mean_serial_interval_real) * Ct[i,t+mean_serial_interval_int] 
          + mean_serial_interval_real * Ct[i,t+mean_serial_interval_int+1]
        );
      }
    }
  }

  { // flux matrices
    fluxt[1,] = to_row_vector(diag_matrix(rep_vector(1.0,N)));
    for (f in 2:F1)
      fluxt[f,] = to_row_vector(flux[f-1]');
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
  // GP prior for transmission rates
  real<lower=0.0> gp_sigma;
  real<lower=0.0,upper=0.632> gp_space_decay;
  real<lower=0.0,upper=0.632> gp_time_decay;
  real<lower=0.0> global_sigma;
  real<lower=0.0> local_space_sigma;
  real<lower=0.0> local_time_sigma;
  vector[N*(Mstep+Mforw)] gp_eta_in;
  vector[(DO_METAPOP && DO_IN_OUT) ? N*(Mstep+Mforw) : 1] gp_eta_out;

  // latent epidemic process
  real<lower=0.0> infection_dispersion;
  vector[N] infection_eta[Tlik];
  real<lower=0.0> xi; // exogeneous infections

  // metapopulation model
  real coupling_mu_eta; // prior mean for sigmoid^-1(coupling_rate) 
  real<lower=0> coupling_sigma; // prior scale for sigmoid^-1(coupling_rate) 
  real<lower=0,upper=1> coupling_alpha1; // 1-autocorrelation for sigmoid^-1(coupling_rate) 
  vector[DO_METAPOP ? (Mstep+Mforw) : 1] coupling_eta; // noise for sigmoid^-1(coupling_rate) 
  simplex[max(1,F)] flux_probs;
  simplex[7] weekly_case_variations;

  // observation process
  vector<lower=0.0>[N] case_precision;
}

transformed parameters {
  // GP prior for transmission rates
  real gp_space_length_scale;
  real gp_time_length_scale;
  matrix[N, Mstep+Mforw] Rin;                 // instantaneous reproduction number
  matrix[N, Mstep+Mforw] Rout;                 // instantaneous reproduction number

  // latent epidemic process
  matrix[N,Tcur] Xt;
  real<lower=0> xi_scale = 0.1;

  // metapopulation model
  row_vector[F1] fluxproportions[Mstep+Mforw];
  real coupling_alpha = 1.0-coupling_alpha1; // autocorrelation for sigmoid^-1(coupling_rate) 
  vector[Mstep+Mforw] coupling_rate;

  { // Prior on transmission rates
    matrix[N,N] K_space;
    matrix[N,N] L_space;
    matrix[Mstep+Mforw,Mstep+Mforw] L_time;
    matrix[Mstep+Mforw,Mstep+Mforw] K_time;

    if (fixed_gp_time_length_scale <= 0.0) {
      gp_time_length_scale = - gp_time_scale / log(1.0-gp_time_decay);
    } else {
      gp_time_length_scale = fixed_gp_time_length_scale;
    }
    K_time  = (
      kernel(TEMPORAL_KERNEL,timedist, 1.0, gp_time_length_scale,
        NONE_KERNEL,EXP_QUAD_KERNEL,MATERN12_KERNEL,MATERN32_KERNEL,MATERN52_KERNEL)
      + diag_matrix(rep_vector(LOCAL_KERNEL ? square(local_time_sigma) : 1e-6, Mstep+Mforw))
    );
    K_time  = K_time .* timecorcut;  // Zero out uncorrelated time entries
    L_time = cholesky_decompose(K_time)';

    if (GLOBAL_KERNEL || SPATIAL_KERNEL != NONE_KERNEL) {
      if (fixed_gp_space_length_scale <= 0.0) {
        gp_space_length_scale = - gp_space_scale / log(1.0-gp_space_decay);
      } else {
        gp_space_length_scale = fixed_gp_space_length_scale;
      }
      K_space = (
        kernel(SPATIAL_KERNEL,geodist, gp_sigma, gp_space_length_scale,
          NONE_KERNEL,EXP_QUAD_KERNEL,MATERN12_KERNEL,MATERN32_KERNEL,MATERN52_KERNEL)
        + diag_matrix(rep_vector(LOCAL_KERNEL ? square(local_space_sigma) : 1e-6, N))
      );
      if (GLOBAL_KERNEL) 
        K_space += rep_matrix(square(global_sigma),N,N);
      L_space = cholesky_decompose(K_space);

      Rin = exp( (L_space * to_matrix(gp_eta_in, N, Mstep+Mforw) * L_time) );
      if (DO_METAPOP && DO_IN_OUT) {
        Rout = exp( (L_space * to_matrix(gp_eta_out, N, Mstep+Mforw) * L_time) );
      } else {
        Rout = rep_matrix(1.0,N,Mstep+Mforw);
      }
 
    } else {
      Rin = exp( local_space_sigma * to_matrix(gp_eta_in, N, Mstep+Mforw) * L_time );
      if (DO_METAPOP && DO_IN_OUT) {
        Rout = exp( (local_space_sigma * to_matrix(gp_eta_out, N, Mstep+Mforw) * L_time) );
      } else {
        Rout = rep_matrix(1.0,N,Mstep+Mforw);
      }
    }

    if (CONSTANT_FORWARD_RT) {
      for (m in (Mstep+1):(Mstep+Mforw)) {
        Rin[:,m] = Rin[,Mstep];
        Rout[:,m] = Rout[,Mstep];
      }
    }
  }


  // metapopulation model
  if (DO_METAPOP) {
    // AR1 process for sigmoid^-1(coupling_rate)
    vector[Mstep+Mforw] Yt;
    real delta = sqrt(1.0-square(coupling_alpha));
    Yt[1] = coupling_eta[1];
    for (i in 2:(Mstep+Mforw))
      Yt[i] = coupling_alpha * Yt[i-1] + delta * coupling_eta[i];
    coupling_rate = inv_logit(
      coupling_mu_loc + 
      coupling_mu_scale * coupling_mu_eta + 
      coupling_sigma * Yt
    );
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

  { // latent renewal process and observation model
    Xt[,1:Tcond] = Xcond;
    for (m in 1:Mstep) {
      matrix[N,N] fluxmatrix;
      if (DO_METAPOP)
        fluxmatrix = to_matrix(fluxproportions[m] * fluxt, N, N);

      for (j in 1:Tstep) {
        int s = (m-1)*Tstep + j;
        int t = Tcond + s;
        int L = min(Tip,t-1);
        vector[N] Zt = Xt[,t-L:t-1] * infprofile_rev[Tip-L+1:Tip];
        vector[N] EXt;
        if (DO_METAPOP) {
          if (DO_IN_OUT)
            EXt = col(Rin,m) .* (xi + fluxmatrix * (Zt .* col(Rout,m)));
          else {
            EXt = col(Rin,m) .* (xi + fluxmatrix * Zt);
          }
        } else {
          EXt = col(Rin,m) .* (xi + Zt);
        }
        Xt[,t] = fabs(
            EXt +
            sqrt((1.0+infection_dispersion) * EXt) .* infection_eta[s]
        );
      }
    }
  }

}

model {
  // prior on transmission rates
  gp_space_decay ~ normal(0.0,gp_space_decay_scale);
  gp_time_decay ~ normal(0.0,gp_time_decay_scale);
  gp_sigma ~ normal(0.0, 0.5);
  global_sigma ~  normal(0.0, 0.5);
  local_space_sigma ~ normal(0.0, 0.5);
  local_time_sigma ~ normal(0.0, 0.1);
  gp_eta_in ~ std_normal();
  gp_eta_out ~ std_normal();
 
  // latent epidemic process 
  infection_dispersion ~ normal(0.0,2.5);
  xi ~ normal(0.0, xi_scale);
  for (s in 1:Tlik)
    infection_eta[s] ~ std_normal();

  // metapopulation model
  flux_probs ~ dirichlet(rep_vector(1.0,max(1,F)));
  coupling_mu_eta ~ std_normal();
  coupling_sigma ~ normal(0.0, coupling_sigma_scale);
  coupling_alpha1 ~ normal(0.0, coupling_alpha_scale);
  coupling_eta ~ std_normal();

  // observation model
  case_precision ~ normal(0.0,10.0);
  weekly_case_variations ~ dirichlet(rep_vector(5.0,7));

  // likelihood
  for (t in Tcond+1:Tcur) {
    int d = (t % 7)+1;
    vector[N] ECt = 
      (7.0*weekly_case_variations[d]) * 
      (Xt[,t-Tdp+1:t] * delayprofile_rev);
    for (j in 1:N) {
      Ct[j,t] ~ neg_binomial_2(ECt[j], ECt[j] * case_precision[j]);
    }
  }

}

generated quantities {

  // Rt
  vector[N] Rt[Mstep+Mproj];
  real Rt_all[Mstep+Mproj];
  // predicted and projected counts
  matrix[N,Tlik] Cpred;
  matrix[N,Tproj] Cproj;
  matrix[N,Tpred] Ppred = rep_matrix(0.0,N,Tpred); // TODO

  { // latent renewal process and observation model
    vector[N] N0 = rep_vector(0.0,N);
    vector[N] N1 = rep_vector(1.0,N);
    matrix[N,Tcur+Tproj] Xproj;

    Xproj[,1:Tcur] = Xt;
    for (m in 1:Mstep+Mproj) {
      matrix[N,N] fluxmatrix;
      vector[N] sum_Zt = rep_vector(0.0,N);
      vector[N] sum_EXt = rep_vector(0.0,N);
      if (DO_METAPOP)
        fluxmatrix = to_matrix(fluxproportions[m] * fluxt, N, N);

      for (j in 1:Tstep) {
        int s = (m-1)*Tstep + j;
        int t = Tcond + s;
        int L = min(Tip,t-1);
        vector[N] EXproj;
        // latent process
        vector[N] Zt;
        Zt = Xproj[,t-L:t-1] * infprofile_rev[Tip-L+1:Tip];
        if (DO_METAPOP) {
          if (DO_IN_OUT)
            EXproj = col(Rin,m) .* (xi + fluxmatrix * (Zt .* col(Rout,m)));
          else {
            EXproj = col(Rin,m) .* (xi + fluxmatrix * Zt);
          }
        } else {
          EXproj = col(Rin,m) .* (xi + Zt);
        }
        if (m > Mstep) {
          Xproj[,t] = fabs(
            EXproj +
            sqrt((1.0+infection_dispersion) * EXproj) .* to_vector(normal_rng(N0,N1))
          );
        }

        // Rt computations
        sum_Zt += Zt;
        sum_EXt += EXproj;

      }
      Rt_all[m] = sum(sum_EXt) / sum(sum_Zt);
      Rt[m] = sum_EXt ./ sum_Zt;
    }

    //observation model
    for (t in Tcond+1:Tcur+Tproj) {
      int s = t - Tcond;
      int d = (t % 7)+1;
      vector[N] ECt = 
        (7.0*weekly_case_variations[d]) *
        Xproj[,t-Tdp+1:t] * delayprofile_rev;
      for (j in 1:N) {
        if (t<=Tcur) {
          Cpred[j,t-Tcond] = ECt[j];
          // Cpred[j,t-Tcond] = neg_binomial_2_rng(ECt[j], case_precision[j]);
        } else {
          Cproj[j,t-Tcur] = ECt[j];
          // Cproj[j,t-Tcur] = neg_binomial_2_rng(ECt[j], case_precision[j]);
        }
      }
    }

  }


  { // print stats
    if (uniform_rng(0.0,1.0)<0.1) {
      print(
        "space ", gp_space_length_scale,
        "; time ", gp_time_length_scale,
        "; sigmas gp ", gp_sigma," g ", global_sigma,
           " l ", local_space_sigma, " ", local_time_sigma,
        "; dispersion X ",infection_dispersion," C ",case_precision[1],
        "; coupling ",coupling_rate[Mstep],
        "; last Rt ", Rt_all[Mstep-1], ", ", Rt_all[Mstep], // ", ", Rt_all[Mstep+Mforw],
        "; last Rin ", Rin[1,Mstep-1], ", ", Rin[1,Mstep] //, ", ", Rin[Mstep+Mforw][1]
      );
    }
  }
}

