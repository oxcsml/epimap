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
      vector Rt, matrix convlik, real coupling_rate, real rad_prob, matrix flux, matrix convrad, matrix convunif) {
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
      vector Rt, matrix convlik, real coupling_rate, real rad_prob, matrix flux, matrix convrad, matrix convunif) {
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
      vector Rt, matrix convlik, real coupling_rate, real rad_prob, matrix flux, matrix convrad, matrix convunif) {
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
      vector Rt, matrix convlik, real coupling_rate, real rad_prob, matrix flux, matrix convrad, matrix convunif) {
    int T = cols(convlik);
    int N = rows(convlik);
    matrix[N,T] convout;
    
    for (i in 1:T) {
      for (j in 1:N) {
        convout[j,i] = Rt[j] * (
            (1.0-coupling_rate) * convlik[j,i] + 
            coupling_rate * convrad[j,i]
        );
      }
    }
    return convout;
  }

  matrix radiation_uniform_in_metapop(
      vector Rt, matrix convlik, real coupling_rate, real rad_prob, matrix flux, matrix convrad, matrix convunif) {
    int T = cols(convlik);
    int N = rows(convlik);
    matrix[N,T] convout;
    
    for (i in 1:T) {
      for (j in 1:N) {
        convout[j,i] = Rt[j] * (
            (1.0-coupling_rate) * convlik[j,i] + 
            coupling_rate * (
                rad_prob * convrad[j,i] + 
                (1.0-rad_prob) * convunif[j,i]
            )
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

  matrix uniform_in_compute_flux(matrix convlik, matrix flux) {
    int T = cols(convlik);
    int N = rows(convlik);
    row_vector[T] convavg;
    matrix[N,T] convflux;

    for (i in 1:T) {
      convavg[i] = 0.0;
      for (j in 1:N)
        convavg[i] += convlik[j,i];
      convavg[i] = convavg[i] / N;
      for (j in 1:N) 
        convflux[j,i] = convavg[i];
    }
    return convflux;
  }

  matrix in_metapop(
      vector Rt, matrix convlik, vector fluxproportions, matrix convflux) {
    int T = cols(convlik);
    int N = rows(convlik);
    int F = cols(convflux);
    // convout[N,T] = diag_pre_multiple(Rt[N],to_matrix(fluxproportions[F] * convflux[F,N*T],N,T ))
    return diag_pre_multiply(Rt,to_matrix(to_row_vector(fluxproportions) * convflux,N,T));

    //matrix[N,T] convout;
    //for (i in 1:T) {
    //  for (j in 1:N) {
    //    vector[F] convin;
    //    for (f in 1:F) {
    //      convin[f] = convflux[f][j,i];
    //    }
    //    convout[j,i] = Rt[j] * dot_product(fluxproportions,convin);
    //  }
    //}
    //return convout;
  }

  matrix compute_flux(matrix convlik, matrix fluxt) {
    int T = cols(convlik);
    int N = rows(convlik);
    int F = rows(fluxt)/N;
    // conflux[F*N,T] = fluxt[F*N,N] * convlik [N,T]
    return to_matrix(fluxt * convlik,F,N*T);

    // matrix[N,T] convflux[F];
    //convflux[1] = convlik;
    //if (F >= 2) {
    //  for (i in 1:T) {
    //    real convin = mean(convlik[,i]);
    //    for (j in 1:N) {
    //      convflux[2][j,i] = convin;
    //    }
    //  }
    //  for (f in 3:F) {
    //    for (i in 1:T) {
    //      for (j in 1:N) {
    //        convflux[f][j,i] = dot_product(flux[f][,j], convlik[,i]);
    //      }
    //    }
    //  }
    //}
    //return convflux;
  }

  // Case count likelihood choices

  real poisson_likelihood_lpmf(int[] count, vector mu, real precision) {
    return poisson_lpmf(count | mu);
  }
  real poisson_predictive_lpmf(int count, real mu, real precision) {
    return poisson_lpmf(count | mu);
  }

  real negative_binomial_2_likelihood_lpmf(int[] count, vector mu, real precision) {
    return neg_binomial_2_lpmf(count | mu, precision);
  }
  real negative_binomial_2_predictive_lpmf(int count, real mu, real precision) {
    return neg_binomial_2_lpmf(count | mu, precision);
  }

  real negative_binomial_3_likelihood_lpmf(int[] count, vector mu, real x) {
    return neg_binomial_2_lpmf(count | mu, mu / x);
  }
  real negative_binomial_3_predictive_lpmf(int count, real mu, real x) {
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
  int<lower=1> F;
  matrix[N,N] flux[F];      // fluxes for metapopulation models. First one for within area infections
}

///////////////////////////////////////////////////////////////////////////

transformed data {
  int Tcur = Tcond+Tlik;    // index of day on which we are estimating Rt
  int Tpred = Tall-Tcur;    // number of days to calculate predictive probabilities for

  vector[D] infprofile_rev; // reversed infection profile

  vector[Tall] Creal[N];     // real type version of Count
  int Clik[N*Tlik];
  int Cpred[N*Tpred];

  matrix[F*N,N] fluxt;

  vector[F] alphas;         // dirichlet parameters for fluxproportions

  // precompute convolutions between Count and infprofile 
  matrix[N,Tlik] convlik;      // for use in likelihood computation
  //matrix[N,Tlik] convlikrad;   // for use in likelihood computation
  //matrix[N,Tlik] convlikunif;  // for use in likelihood computation
  matrix[F,N*Tlik] convlikflux;
  matrix[N,Tpred] convpred;    // for use in predictive probs of future counts
  matrix[F,N*Tpred] convpredflux;
  //matrix[N,Tpred] convpredrad; // for use in predictive probs of future counts
  //matrix[N,Tpred] convpredunif;// for use in predictive probs of future counts
  matrix[N,Tproj] convproj;    // for use in forecasting into future 

  // reverse infection profile
  for (i in 1:D)
    infprofile_rev[i] = infprofile[D-i+1];

  {
    matrix[F,N*N] fluxtmp;
    for (f in 1:F) 
      fluxtmp[f,] = to_row_vector(flux[f]');
    fluxt = to_matrix(fluxtmp,F*N,N);
  }

  alphas[1] = F;
  for (f in 2:F)
    alphas[f] = 1.0;

  Clik = to_array_1d(Count[,Tcond+1:Tcond+Tlik]);
  Cpred = to_array_1d(Count[,Tcond+Tlik+1:Tcond+Tlik+Tpred]);

  for (j in 1:N) {
    Creal[j] = to_vector(Count[j,]);
  //  for (i in 1:Tall)
  //    Creal[j,i] = Count[j,i];

    // precompute convolutions between counts and infprofile
    // 1e-6 for numerical issues when case counts 0
    for (i in 1:Tlik) {
      int L = min(D,Tcond+i-1); // length of infection profile that overlaps with case counts 
      convlik[j,i] = dot_product(Creal[j][Tcond-L+i:Tcond-1+i], infprofile_rev[D-L+1:D]) + 1e-6;
    }
    for (i in 1:Tpred) {
      int L = min(D,Tcur+i-1); // length of infection profile that overlaps with case counts 
      convpred[j,i] = dot_product(Creal[j][Tcur-L+i:Tcur-1+i], infprofile_rev[D-L+1:D]) + 1e-6;
    }
    for (i in 1:Tproj) {
      int L = min(D,Tcur+i-1); // length of infection profile that overlaps with case counts 
      convproj[j,i] = dot_product(Creal[j][Tcur-L+i:Tcur], infprofile_rev[D-L+1:D-i+1]) + 1e-6;
    }
  }
  convlikflux = compute_flux(convlik,fluxt);
  convpredflux = compute_flux(convpred,fluxt);
  //convlikunif = uniform_in_compute_flux(convlik,flux);
  //convpredunif = uniform_in_compute_flux(convpred,flux);
  //convlikrad = radiation_in_compute_flux(convlik,flux);
  //convpredrad = radiation_in_compute_flux(convpred,flux);
}

///////////////////////////////////////////////////////////////////////////

parameters {
  real<lower=0> gp_length_scale;
  real<lower=0> gp_sigma;
  real<lower=0> global_sigma;
  real<lower=0> local_exp[N];
  real<lower=0> local_scale;
  vector[N] eta;

  real<lower=0> precision;
  simplex[F] fluxproportions;
}


///////////////////////////////////////////////////////////////////////////

transformed parameters {
  vector[N] Rt;                 // instantaneous reproduction number
  real<lower=0> local_sigma2[N];

  {
    matrix[N,N] K;
    matrix[N,N] L;
    real global_sigma2 = square(global_sigma);

    for (i in 1:N) 
      local_sigma2[i] = (2.0 * square(local_scale)) * local_exp[i];
 
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
  vector[N*Tlik] convout;

  // Ravg ~ normal(1.0,1.0);
  precision ~ normal(0.0,10.0);

  // GP prior density
  eta ~ std_normal();
  // gp_length_scale ~ normal(0.0,5.0);
  gp_length_scale ~ gig(2, 2.0, 2.0);
  gp_sigma ~ normal(0.0, 0.25);
  global_sigma ~ normal(0.0, 0.25);
  local_scale ~ normal(0.0, 0.1);
  for (i in 1:N) {
    local_exp[i] ~ exponential(1.0);
    // local_sigma2[i] ~ exponential(0.5 / square(local_scale)); // reparameterised
  }
  fluxproportions ~ dirichlet(alphas);

  // metapopulation infection rate model
  // convout = METAPOP_metapop(Rt,convlik,coupling_rate,rad_prob,flux,convlikrad,convlikunif);
  convout = to_vector(in_metapop(Rt,convlik,fluxproportions,convlikflux));

  // compute likelihoods
  Clik ~ OBSERVATION_likelihood(convout, precision);
  //for (j in 1:N) {
  //  for (i in 1:Tlik) {
  //    Clik[j,i] ~ OBSERVATION_predictive(convout[j,i], precision);
  //  }
  //}

}



///////////////////////////////////////////////////////////////////////////

generated quantities {
  real R0;
  matrix[N,Tpred] Ppred;
  matrix[N,Tproj] Cproj; 

  // Estimated R0 over all areas
  {
    //real convsum = sum(convlik);
    R0 = sum(to_row_vector(Rt) * convlik) / sum(convlik);
    //for (i in 1:Tlik) {
    //  for (j in 1:N) {
    //    R0 += Rt[j] * convlik[j,i];
    //    convsum += convlik[j,i];
    //  }
    //}
    //R0 = R0 / convsum;
  }

  // predictive probability of future counts
  {
    //matrix[N,Tpred] convout = METAPOP_metapop(Rt,convpred,coupling_rate,rad_prob,flux,convpredrad,convpredunif);
    matrix[N,Tpred] convout = in_metapop(Rt,convpred,fluxproportions,convpredflux);
    for (i in 1:Tpred)
      for (j in 1:N)
        Ppred[j,i] = exp(OBSERVATION_predictive_lpmf(Count[j,Tcur+i] |
            convout[j,i], precision));
  }

  // forecasting *mean* counts given parameters
  {
    matrix[N,1] convprojall;
    //matrix[N,1] convprojrad;
    //matrix[N,1] convprojunif;
    matrix[F,N*1] convprojflux;
    for (i in 1:Tproj) {
      for (j in 1:N) 
        convprojall[j,1] = convproj[j,i] + 
            dot_product(Cproj[j][1:(i-1)], infprofile_rev[(D-i+2):D]);
      convprojflux = compute_flux(convprojall,fluxt);
      //convprojrad  = radiation_in_compute_flux(convprojall,flux);
      //convprojunif = uniform_in_compute_flux(convprojall,flux);
      //Cproj[:,i] = METAPOP_metapop(Rt,convprojall,coupling_rate,rad_prob,flux,convprojrad,convprojunif)[:,1];
      Cproj[:,i] = in_metapop(Rt,convprojall,fluxproportions,convprojflux)[:,1];
    }
  }

}


