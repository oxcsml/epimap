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
  int<lower=1> Tall;        // number of all days in case count time series
  int<lower=1> Tcond;       // number of days we will condition on
  int<lower=1> Tstep;       // number of days to step for each time step of Rt prediction
  int<lower=1> Tpred;
  int<lower=1> Tproj;       // number of days to forecast
  int Count[N, Tall];       // case counts,

  int<lower=1> Tip;         // length of infection profile
  vector[Tip] infprofile;   // infection profile aka serial interval distribution
  int<lower=1> Tdp;
  vector[Tdp] delayprofile; // delay profile

  matrix[N,N] geodist;      // distance between locations
  matrix[M,M] timedist;     // distance between time samples
  matrix[M,M] timecorcut;   // matrix specifying which time points should be correlated (to account for lockdown)
  int<lower=0,upper=1> do_metapop;
  int<lower=0,upper=1> do_in_out;
  int F;
  matrix[N,N] flux[F];      // fluxes for radiation metapopulation model
}

transformed data {
  int Tlik = M*Tstep;
  int Tcur = Tall-Tpred;    // index of last day modelled

  matrix[N,Tcond] Ccond;     // condtioned counts
  vector[Tip] infprofile_rev; // reversed infection profile
  vector[Tdp] delayprofile_rev; // reversed infection profile

  int F1 = F+1;
  vector[max(1,F)] ones = rep_vector(1.0,max(1,F));
  matrix[F1,N*N] fluxt;      // transposed flux matrices

  // reverse infection and delay profile
  for (i in 1:Tip)
    infprofile_rev[i] = infprofile[Tip-i+1];
  for (i in 1:Tdp)
    delayprofile_rev[i] = delayprofile[Tdp-i+1];

  Ccond = rep_matrix(1e-2,N,Tcond);
  for (i in 1:N) 
    for (t in 1:Tcond) 
      Ccond[i,t] += Count[i,t];

  {
    fluxt[1,] = to_row_vector(diag_matrix(rep_vector(1.0,N)));
    for (f in 2:F1)
      fluxt[f,] = to_row_vector(flux[f-1]');
  }

}

parameters {
  real<lower=0> gp_time_length_scale;
  real<lower=0> gp_space_length_scale;
  real<lower=0> gp_space_sigma;

    
  matrix<lower=0>[N,M] local_exp;
  real<lower=0> local_scale;
  real<lower=0> global_sigma;

  vector[N*M] eta_in;
  vector[N*M] eta_out;

  vector[N*M] epsilon_in;
  vector[N*M] epsilon_out;

  vector[N] Ceta[Tlik];

  real<lower=0> phi_latent;
  real<lower=0> phi_observed;
  real<lower=0,upper=1> coupling_rate[M];
  simplex[max(1,F)] flux_probs;
}

transformed parameters {
  matrix[N, M] Rin;                 // instantaneous reproduction number
  matrix[N, M] Rout;                 // instantaneous reproduction number
  matrix[N, M] Rt;
  vector[M]    R0;

  matrix[N, M] local_sigma;

  row_vector[F1] fluxproportions[M];

  matrix[N,Tlik] Cinfer;
  matrix[N,Tlik] Cpred;
  
  { // GP
    matrix[N,N] K_space;
    matrix[M,M] K_time;
    matrix[N,N] L_space;
    matrix[M,M] L_time;
    real global_sigma2 = square(global_sigma);

    // GP space kernel
    K_space = SPATIAL_kernel(geodist, gp_space_sigma, gp_space_length_scale); // space kernel
    K_space += GLOBAL_var(global_sigma2); // Add global noise

    // GP time kernel
    K_time = TEMPORAL_kernel(timedist, 1.0, gp_time_length_scale);
    K_time .*= timecorcut;  // Zero out uncorrelated time entries

    L_space = cholesky_decompose(K_space);
    L_time = cholesky_decompose(K_time);

    // Noise kernel (reshaped from (N*M, N*M) to (N,M) since diagonal)
    for (m in 1:M) {
      for (n in 1:N) {
        local_sigma[n, m] = sqrt(LOCAL_var(local_exp[n, m] * (2 * square(local_scale))));
      }
    }


    // Rt 
    // Compute (K_time (*) K_space) * eta via efficient kronecker trick. Don't reshape back to vector for convinience.
    // Add on the location-time dependant noise as well
    Rin = exp(
        (L_space * to_matrix(eta_in, N, M) * L_time') + 
        (local_sigma .* to_matrix(epsilon_in, N, M))
    );
    if (do_metapop && do_in_out) {
      Rout = exp(
          (L_space * to_matrix(eta_out, N, M) * L_time') + 
          (local_sigma .* to_matrix(epsilon_out, N, M))
      );
    } else {
      Rout = rep_matrix(1.0,N,M);
    }
  }

  { // metapopulation infection rate model
    if (do_metapop) {
      for (m in 1:M){
        fluxproportions[m, 1] = 1.0-coupling_rate[m];
      for (f in 2:F1)
        fluxproportions[m, f] = coupling_rate[m]*flux_probs[f-1];
      }
    } else {
      for (m in 1:M){
        fluxproportions[m, 1] = 1.0;
        for (f in 2:F1)
          fluxproportions[m, f] = 0.0;
      }
      
    }
  }

  { // latent renewal process and observation model
    matrix[N,Tcur] Clatent;
    Clatent[,1:Tcond] = Ccond;
    for (m in 1:M) {
      matrix[N,N] fluxmatrix;
      vector[N] sum_local_infforce = rep_vector(0.0,N);
      vector[N] sum_total_infforce = rep_vector(0.0,N);
      if (do_metapop) 
        fluxmatrix = to_matrix(fluxproportions[m] * fluxt, N, N);

      for (j in 1:Tstep) {
        int s = (m-1)*Tstep + j;
        int t = Tcond + s;
        int L = min(Tip,t-1);
        // latent process
        vector[N] local_infforce;
        vector[N] total_infforce;
        local_infforce = Clatent[,t-L:t-1] * infprofile_rev[Tip-L+1:Tip];
        total_infforce = local_infforce;
        if (do_metapop) {
          if (do_in_out) 
            total_infforce .*= col(Rout,m); 
          total_infforce = fluxmatrix * total_infforce;
        }
        total_infforce .*= col(Rin,m);
        Clatent[,t] = 1e-2 + fabs(
            total_infforce + 
            sqrt((1.0+phi_latent) * total_infforce) .* Ceta[s]
        );

        // Rt computations
        sum_local_infforce += local_infforce;
        sum_total_infforce += total_infforce;

        // observation model
        Cpred[,s] = Clatent[,t-Tdp:t-1] * delayprofile_rev;
      }
      R0[m] = sum(sum_total_infforce) / sum(sum_local_infforce);
      Rt[,m] = sum_total_infforce ./ sum_local_infforce;
    }
    Cinfer = Clatent[,Tcond+1:Tcur];
  }

}

model {
  coupling_rate ~ normal(0.0, .1);
  flux_probs ~ dirichlet(ones);
  phi_latent ~ normal(0.0,1.0);
  phi_observed ~ normal(0.0,5.0);

  // GP prior density
  eta_in ~ std_normal();
  eta_out ~ std_normal();
  epsilon_in ~ std_normal();
  epsilon_out ~ std_normal();

  gp_time_length_scale ~ gig(7, 0.2, 1.0);
  gp_space_length_scale ~ gig(5, 5.0, 5.0);
  gp_space_sigma ~ normal(0.0, 0.25);
  global_sigma ~  normal(0.0, 0.25);
  local_scale ~ normal(0.0, 0.1);
  for (j in 1:M){
    for (i in 1:N) {
      local_exp[i, j] ~ exponential(1.0);
    }
  }

  //print("phi_observed: ",phi_observed);
  //print("phi_latent: ",phi_latent);
  //print("coupling_rate: ",coupling_rate);
  //print("gp_space_sigma: ",gp_space_sigma);
  //print("gp_space_length_scale: ",gp_space_length_scale);
  //print("gp_time_length_scale: ",gp_time_length_scale);
  //print("global_sigma: ",global_sigma);
  //print("local_scale: ",local_scale);

  for (s in 1:Tlik) 
    Ceta[s] ~ std_normal();

  { // compute likelihoods
    real err = 0.0;
    for (t in Tcond+1:Tcur) {
      int s = t - Tcond;
      for (j in 1:N) {
        //Count[j,t] ~ OBSERVATION_likelihood(Ecase[j], phi_observed);
        Count[j,t] ~ neg_binomial_2(Cpred[j,s], Cpred[j,s] / (1e-6+phi_observed));
      }
    }
  }
}

generated quantities {
  real err = 0.0;
  {
    for (t in Tcond+1:Tcur) {
      int s = t - Tcond;
      for (j in 1:N) {
        err += fabs(Cpred[j,s] - Count[j,t]);
      }
    }
    err = err/Tlik/N;
  }

  if (uniform_rng(0.0,1.0)<.1) {
    int ind[30];
    for (i in 1:30) ind[i] = 10*i;
    print("Clatent: ",round(Cinfer[ind,Tlik]));
    print("phi_observed: ",phi_observed);
    print("error: ",err);
  }

}




