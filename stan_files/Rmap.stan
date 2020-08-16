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

 
  matrix none_metapop(vector Rin, matrix convlik) {
    return diag_pre_multiply(Rin, convlik);
  }

  matrix in_compute_flux(matrix convlik, matrix fluxt) {
    int T = cols(convlik);
    int N = rows(convlik);
    int F1 = rows(fluxt);
    // convflux[F1,N*T] = to_matrix(to_matrix(fluxt[F1,N*N], F1*N, N) * convlik[N,T], F1,N*T)
    return to_matrix(to_matrix(fluxt, F1*N, N) * convlik, F1, N*T);
  }

  matrix in_metapop(
      vector Rin, matrix convlik, 
      row_vector fluxproportions, matrix convflux) {
    int T = cols(convlik);
    int N = rows(convlik);
    // convout[N,T] = diag_pre_multiply(Rin[N],
    //     to_matrix(fluxproportions[F1] * convflux[F1,N*T], N, T))
    return diag_pre_multiply(Rin,
        to_matrix(fluxproportions * convflux, N, T));
  }

  matrix in_out_metapop(
      vector Rin, vector Rout, matrix convlik,
      row_vector fluxproportions, matrix fluxt) {
    int N = rows(convlik);
    // convin[N,T] = diag_pre_multiply(Rin[N], 
    //     to_matrix(fluxproportions[F] * fluxt[F,N*N], N, N) *
    //     diag_pre_multiply(Rout[N],convlik[N,T])
    // );
    return diag_pre_multiply(Rin, 
        to_matrix(fluxproportions * fluxt, N, N) *
        diag_pre_multiply(Rout,convlik)
    );
  }

  matrix[] in_out_timevary_metapop(
      vector[] Rin, vector[] Rout, matrix[] convlik, row_vector[] fluxproportions, matrix fluxt) {
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
      convin[m] = diag_pre_multiply(Rin[m], 
        to_matrix(fluxproportions[m] * fluxt, N, N) *
        diag_pre_multiply(Rout[m],convlik[m]));
    }
    return convin;
  }

  matrix metapop(
      int do_metapop, int do_in_out,
      vector Rin, vector Rout, matrix convlik, matrix convflux,
      row_vector fluxproportions, matrix fluxt) {
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
  int<lower=0,upper=1> do_metapop;
  int<lower=0,upper=1> do_in_out;
  int F;
  matrix[N,N] flux[F];         // fluxes for radiation metapopulation model
  //matrix[N,N] metaflux[F];
}

///////////////////////////////////////////////////////////////////////////

transformed data {
  int Tcur = Tcond+Tlik;    // index of day on which we are estimating Rt
  int Tpred = Tall-Tcur;    // number of days to calculate predictive probabilities for
  int F1 = F+1;
  vector[max(1,F)] ones = rep_vector(1.0,max(1,F));

  vector[Tall] Creal[N];     // real type version of Count
  vector[D] infprofile_rev; // reversed infection profile

  // precompute convolutions between Count and infprofile 
  matrix[N,Tlik] convlik;      // for use in likelihood computation
  matrix[N,Tlik] convlikrad;   // for use in likelihood computation
  matrix[N,Tlik] convlikunif;  // for use in likelihood computation
  matrix[N,Tpred] convpred;    // for use in predictive probs of future counts
  matrix[N,Tpred] convpredrad; // for use in predictive probs of future counts
  matrix[N,Tpred] convpredunif;// for use in predictive probs of future counts
  matrix[N,Tproj] convproj;    // for use in forecasting into future 

  matrix[F1,N*N] fluxt;      // transposed flux matrices

  matrix[F1,N*Tlik] convlikflux;
  matrix[F1,N*Tpred] convpredflux;
  int Clik[Tlik,N];

  // reverse infection profile
  for (i in 1:D)
    infprofile_rev[i] = infprofile[D-i+1];

  {
    fluxt[1,] = to_row_vector(diag_matrix(rep_vector(1.0,N)));
    for (f in 2:F1)
      fluxt[f,] = to_row_vector(flux[f-1]');
  }

  for (j in 1:N)
    for (i in 1:Tlik)
      Clik[i,j] = Count[j,Tcond+i]; // transposed, vectorises correctly

  for (j in 1:N) {
    for (i in 1:Tall)
      Creal[j,i] = Count[j,i];

    // precompute convolutions between counts and infprofile
    for (i in 1:Tlik) {
      int L = min(D,Tcond+i-1); // length of infection profile that overlaps with case counts 
      convlik[j,i] = dot_product(Creal[j][Tcond-L+i:Tcond-1+i], infprofile_rev[D-L+1:D])+1e-6;
    }
    for (i in 1:Tpred) {
      int L = min(D,Tcur+i-1); // length of infection profile that overlaps with case counts 
      convpred[j,i] = dot_product(Creal[j][Tcur-L+i:Tcur-1+i], infprofile_rev[D-L+1:D])+1e-6;
    }
    for (i in 1:Tproj) {
      int L = min(D,Tcur+i-1); // length of infection profile that overlaps with case counts 
      convproj[j,i] = dot_product(Creal[j][Tcur-L+i:Tcur], infprofile_rev[D-L+1:D-i+1])+1e-6;
    }
  }

  if (do_metapop && !do_in_out) {
    convlikflux = in_compute_flux(convlik,fluxt);
    convpredflux = in_compute_flux(convpred,fluxt);
  }

}

///////////////////////////////////////////////////////////////////////////

parameters {
  real<lower=0> gp_length_scale;
  real<lower=0> gp_sigma;
  real<lower=0> global_sigma;
  real<lower=0> local_exp[N];
  real<lower=0> local_scale;
  vector[N] eta_in;
  vector[N] eta_out;

  real<lower=0> precision;
  // real<lower=0> Ravg;
  real<lower=0,upper=1> coupling_rate;
  simplex[max(1,F)] flux_probs;
}


///////////////////////////////////////////////////////////////////////////

transformed parameters {
  vector[N] Rin;                 // instantaneous reproduction number
  vector[N] Rout;                 // instantaneous reproduction number
  real<lower=0> local_sigma2[N];
  row_vector[F1] fluxproportions;
  matrix[N,Tlik] convlikout;

  {
    matrix[N,N] K;
    matrix[N,N] L;
    real global_sigma2 = square(global_sigma);

    for (i in 1:N) {
      local_sigma2[i] = local_exp[i] * (2.0 * square(local_scale));
    }
 
    // GP prior.
    K = SPATIAL_kernel(geodist, gp_sigma, gp_length_scale); // kernel
    for (i in 1:N) {
      K[i,i] = K[i,i] + LOCAL_var(local_sigma2[i]);
    }
    K = K + GLOBAL_var(global_sigma2);

    L = cholesky_decompose(K);
    Rin = exp(L * eta_in);
    Rout = exp(L * eta_out);
  }

  // metapopulation infection rate model
  if (do_metapop) {
    fluxproportions[1] = 1.0-coupling_rate;
    for (f in 2:F1)
      fluxproportions[f] = coupling_rate*(1.0-flux_probs[f-1]);
  } else {
    fluxproportions[1] = 1.0;
    for (f in 2:F1)
      fluxproportions[f] = 0.0;
  }
  convlikout = metapop(do_metapop,do_in_out,
      Rin,Rout,convlik,convlikflux,fluxproportions,fluxt);
}


///////////////////////////////////////////////////////////////////////////

model {
  // Ravg ~ normal(1.0,1.0);
  coupling_rate ~ normal(0.0, .5);
  flux_probs ~ dirichlet(ones);
  precision ~ normal(0.0,10.0);

  // GP prior density
  eta_in ~ std_normal();
  eta_out ~ std_normal();
  // gp_length_scale ~ normal(0.0,5.0);
  gp_length_scale ~ gig(2, 2.0, 2.0);
  gp_sigma ~ normal(0.0, 0.5);
  global_sigma ~ normal(0.0, 0.5);
  local_scale ~ normal(0.0, 0.5);
  for (i in 1:N) {
    local_exp[i] ~ exponential(1.0);
    // local_sigma2[i] ~ exponential(0.5 / square(local_scale)); // reparameterised
  }

  // compute likelihoods
  for (i in 1:Tlik) {
    for (j in 1:N) {
      Count[j,Tcond+i] ~ OBSERVATION_likelihood(convlikout[j,i], precision);
      //Clik[i,j] ~ OBSERVATION_likelihood(convlikout[j,i], precision);
    }
  }

}



///////////////////////////////////////////////////////////////////////////

generated quantities {
  real R0;
  vector[N] Rt;
  matrix[N,Tpred] Ppred;
  matrix[N,Tproj] Cproj; 

  // Estimated R0 and Rt for all areas
  {
    vector[N] oneN = rep_vector(1.0,N);
    vector[Tlik] oneT = rep_vector(1.0,Tlik);
    matrix[N,Tlik] convone = metapop(do_metapop,do_in_out,
        oneN,oneN,convlik,convlikflux,fluxproportions,fluxt);
    R0 = sum(convlikout) / sum(convone);
    Rt = (convlikout * oneT) ./ (convone * oneT);
  }

  // predictive probability of future counts
  {
    //matrix[N,Tpred] convpredout = METAPOP_metapop(Rin,convpred,coupling_rate,rad_prob,flux,convpredrad,convpredunif);
    //matrix[N,Tpred] convpredout = in_metapop(Rin,convpred,fluxproportions,convpredflux);
    matrix[N,Tpred] convpredout = metapop(do_metapop,do_in_out,
        Rin,Rout,convpred,convpredflux,fluxproportions,fluxt);
    for (i in 1:Tpred)
      for (j in 1:N)
        Ppred[j,i] = exp(OBSERVATION_likelihood_lpmf(Count[j,Tcur+i] |
            convpredout[j,i], precision));
  }

  // forecasting *mean* counts given parameters
  {
    matrix[N,1] convprojall;
    //matrix[N,1] convprojrad;
    //matrix[N,1] convprojunif;
    matrix[F1,N*1] convprojflux;
    for (i in 1:Tproj) {
      for (j in 1:N) 
        convprojall[j,1] = convproj[j,i] + 
            dot_product(Cproj[j][1:(i-1)], infprofile_rev[(D-i+2):D]);
      if (do_metapop && !do_in_out)
        convprojflux = in_compute_flux(convprojall,fluxt);
      Cproj[:,i] = metapop(do_metapop,do_in_out,
          Rin,Rout,convprojall,convprojflux,fluxproportions,fluxt)[:,1];
    }
  }

}


