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
      vector Rout, matrix convlik, real coupling_rate, real rad_prob, matrix flux, matrix convrad, matrix convunif) {
    int T = cols(convlik);
    int N = rows(convlik);
    row_vector[T] convavg;
    matrix[N,T] convout;
    
    for (i in 1:T) {
      convavg[i] = 0.0;
      for (j in 1:N)
        convavg[i] += Rout[j] * convlik[j,i];
      convavg[i] = convavg[i] / N;
    }

    for (j in 1:N)
      for (i in 1:T)
        convout[j,i] = (
            (1.0-coupling_rate) *  Rout[j] * convlik[j,i] +
            coupling_rate * convavg[i]
        );
    return convout;
  }

  matrix uniform_in_metapop(
      vector Rin, matrix convlik, real coupling_rate, real rad_prob, matrix flux, matrix convrad, matrix convunif) {
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
        convout[j,i] = Rin[j] * (
            (1.0-coupling_rate) * convlik[j,i] +
            coupling_rate * convavg[i]
        );
    return convout;
  }

  matrix radiation_out_metapop(
      vector Rout, matrix convlik, real coupling_rate, real rad_prob, matrix flux, matrix convrad, matrix convunif) {
    int T = cols(convlik);
    int N = rows(convlik);
    matrix[N,T] convout;
    
    for (i in 1:T) {
      for (j in 1:N) {
        convout[j,i] = Rout[j] * (1.0-coupling_rate) * convlik[j,i];
        for (k in 1:N) {
          convout[j,i] += Rout[k] * coupling_rate * convlik[k,i] * flux[k,j];
        } 
      }
    }
    return convout;
  }

  matrix radiation_in_metapop(
      vector Rin, matrix convlik, real coupling_rate, real rad_prob, matrix flux, matrix convrad, matrix convunif) {
    int T = cols(convlik);
    int N = rows(convlik);
    matrix[N,T] convout;
    
    for (i in 1:T) {
      for (j in 1:N) {
        convout[j,i] = Rin[j] * (
            (1.0-coupling_rate) * convlik[j,i] + 
            coupling_rate * convrad[j,i]
        );
      }
    }
    return convout;
  }

  matrix radiation_uniform_in_metapop(
      vector Rin, matrix convlik, real coupling_rate, real rad_prob, matrix flux, matrix convrad, matrix convunif) {
    int T = cols(convlik);
    int N = rows(convlik);
    matrix[N,T] convout;
    
    for (i in 1:T) {
      for (j in 1:N) {
        convout[j,i] = Rin[j] * (
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

  matrix new_metapop(
      vector Rin, matrix convlik, row_vector fluxproportions, matrix convflux) {
    int T = cols(convlik);
    int N = rows(convlik);
    int F1 = cols(convflux);
    // convout[N,T] = diag_pre_multiply(Rin[N],
    //     to_matrix(fluxproportions[F1] * convflux[F1,N*T], N, T))
    return diag_pre_multiply(Rin,
        to_matrix(fluxproportions * convflux, N, T));
  }



  matrix none_metapop(
      vector Re, matrix convlik, real coupling_rate, real rad_prob, matrix flux, matrix convrad, matrix convunif) {
    int T = cols(convlik);
    int N = rows(convlik);
    matrix[N,T] convout;

    for (j in 1:N)
      for (i in 1:T)
        convout[j,i] = Re[j] * convlik[j,i];

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


  matrix radiation_in_compute_flux(matrix convlik, matrix flux) {
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


  matrix new_compute_flux(matrix convlik, matrix fluxt) {
    int T = cols(convlik);
    int N = rows(convlik);
    int F1 = rows(fluxt)/N;
    // convflux[F1,N*T] = to_matrix(fluxt[F1*N,N] * convlik[N,T], F1,N*T)
    return to_matrix(fluxt * convlik, F1, N*T);
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
  int F;
  matrix[N,N] flux[F];         // fluxes for radiation metapopulation model
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

  matrix[F1*N,N] fluxt;      // transposed flux matrices

  matrix[F1,N*Tlik] convflux;
  matrix[F1,N*Tpred] convpredflux;
  int Clik[Tlik,N];

  // reverse infection profile
  for (i in 1:D)
    infprofile_rev[i] = infprofile[D-i+1];

  {
    matrix[F1,N*N] fluxtmp;
    fluxtmp[1,] = to_row_vector(diag_matrix(rep_vector(1.0,N)));
    for (f in 2:F1)
      fluxtmp[f,] = to_row_vector(flux[f-1]');
    fluxt = to_matrix(fluxtmp,F1*N,N);
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
  //convlikunif = uniform_in_compute_flux(convlik,flux);
  //convpredunif = uniform_in_compute_flux(convpred,flux);
  //convlikrad = radiation_in_compute_flux(convlik,flux);
  //convpredrad = radiation_in_compute_flux(convpred,flux);

  convflux = new_compute_flux(convlik,fluxt);
  convpredflux = new_compute_flux(convpred,fluxt);
  //{
  //  matrix[N,Tlik] radflux = to_matrix(convflux[3,],N,Tlik);
  //  matrix[N,Tlik] unifflux = to_matrix(convflux[2,],N,Tlik);
  //  matrix[N,Tlik] localflux = to_matrix(convflux[1,],N,Tlik);
  //  print("Check for flux computations: ",
  //        max(localflux-convlik),"  ",min(localflux-convlik),"  ",
  //        max(unifflux-convlikunif),"  ",min(unifflux-convlikunif),"  ",
  //        max(radflux-convlikrad),"  ",min(radflux-convlikrad));
  //  for (j in 1:20) {
  //    print(convlikunif[j,1],"  ",unifflux[j,1]);
  //    print(convlikunif[j,3],"  ",unifflux[j,3]);
  //    print(convlik[j,3],"  ",localflux[j,3]);
  //  }
  //  for (i in 1:Tlik)
  //    print(convlikunif[100,i],"  ",unifflux[100,i]);
  //}
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
  // real<lower=0> Ravg;
  real<lower=0,upper=1> coupling_rate;
  simplex[max(1,F)] flux_probs;
}


///////////////////////////////////////////////////////////////////////////

transformed parameters {
  vector[N] Rin;                 // instantaneous reproduction number
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
    Rin = exp(L * eta);

    fluxproportions[1] = 1.0-coupling_rate;
    for (f in 2:F1) {
      fluxproportions[f] = coupling_rate*(1.0-flux_probs[f-1]);
    }
  }

  convlikout = new_metapop(Rin,convlik,fluxproportions,convflux);

}


///////////////////////////////////////////////////////////////////////////

model {
  vector[Tlik] coupling;

  // Ravg ~ normal(1.0,1.0);
  coupling_rate ~ normal(0.0, .5);
  flux_probs ~ dirichlet(ones);
  precision ~ normal(0.0,10.0);

  // GP prior density
  eta ~ std_normal();
  // gp_length_scale ~ normal(0.0,5.0);
  gp_length_scale ~ gig(2, 2.0, 2.0);
  gp_sigma ~ normal(0.0, 0.5);
  global_sigma ~ normal(0.0, 0.5);
  local_scale ~ normal(0.0, 0.5);
  for (i in 1:N) {
    local_exp[i] ~ exponential(1.0);
    // local_sigma2[i] ~ exponential(0.5 / square(local_scale)); // reparameterised
  }
  // metapopulation infection rate model

  // both below computes same thing.
  //convout = METAPOP_metapop(Rin,convlik,coupling_rate,rad_prob,flux,convlikrad,convlikunif);

  //{
    //matrix[N,Tlik] oldconvout = METAPOP_metapop(Rin,convlik,coupling_rate,rad_prob,flux,convlikrad,convlikunif);
    //print("Check for metapop computations: ",
    //      max(convout-oldconvout),"  ",min(convout-oldconvout));
    //for (i in 1:7) {
    //  print(convout[i,2],"  ",newconvout[i,2]);
    //  print(convout[200,i],"  ",newconvout[200,i]);
    //}
  //}


  // compute likelihoods
  for (i in 1:Tlik) {
    for (j in 1:N) {
      Count[j,Tcond+i] ~ OBSERVATION_likelihood(convlikout[j,i], precision);
      //Clik[i,j] ~ OBSERVATION_likelihood(convout[j,i], precision);
      //print(OBSERVATION_likelihood_lpmf(Count[j,Tcond+i]|convout[j,i], precision),
      //      OBSERVATION_likelihood_lpmf(Clik[i,j]|convout[j,i], precision));

    }
  }

}



///////////////////////////////////////////////////////////////////////////

generated quantities {
  real R0;
  vector[N] Rt;
  matrix[N,Tpred] Ppred;
  matrix[N,Tproj] Cproj; 

  // Estimated R0 over all areas
  {
    vector[N] R1 = rep_vector(1.0,N);
    vector[Tlik] T1 = rep_vector(1.0,Tlik);
    matrix[N,Tlik] convone = new_metapop(R1,convlik,fluxproportions,convflux);
    R0 = sum(convlikout) / sum(convone);
    Rt = (convlikout * T1) ./ (convone * T1);
  }

  // predictive probability of future counts
  {
    //matrix[N,Tpred] convpredout = METAPOP_metapop(Rin,convpred,coupling_rate,rad_prob,flux,convpredrad,convpredunif);
    matrix[N,Tpred] convpredout = new_metapop(Rin,convpred,fluxproportions,convpredflux);
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
      //convprojrad  = radiation_in_compute_flux(convprojall,flux);
      //convprojunif = uniform_in_compute_flux(convprojall,flux);
      convprojflux = new_compute_flux(convprojall,fluxt);
      //Cproj[:,i] = METAPOP_metapop(Rin,convprojall,coupling_rate,rad_prob,flux,convprojrad,convprojunif)[:,1];
      Cproj[:,i] = new_metapop(Rin,convprojall,fluxproportions,convprojflux)[:,1];
    }
  }

}


