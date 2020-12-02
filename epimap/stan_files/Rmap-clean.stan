data {
  int<lower=1> Tall;  //
  int<lower=1> Tcond;  //
  int<lower=1> Tstep; // length of step (7 days = 1 week)
  int<lower=0> Tpred; // number of days to compute prediction likelihood for
  int<lower=1> Nstep; // number of steps
  int<lower=0> Nproj;
  int Count[1,Tall];
  int<lower=1> Tip;
  vector[Tip] infprofile;
  int<lower=1> Ttdp;
  vector[Ttdp] testdelayprofile;
  int<lower=1> Trdp;
  vector<lower=0>[Trdp+1] resultdelayalpha;
  real<lower=0> gp_time_scale;
  real<lower=0> gp_time_decay_scale;
  real fixed_gp_time_length_scale;
  real<lower=0> mu_scale;
  real<lower=0> sigma_scale;
  real<lower=0> phi_latent_scale;
  real<lower=0> phi_observed_scale;
  real<lower=0,upper=1> outlier_prob_threshold;
  int<lower=1> outlier_count_threshold;
  real<lower=0> xi_scale;
  int<lower=0,upper=1> reconstruct_infections;
}

transformed data {
  int Tlik = Nstep*Tstep;
  int Tcur = Tcond+Tlik;
  vector[Tip] infprofile_rev;
  vector[Ttdp] testdelayprofile_rev;
  int mean_serial_interval_int = 0;
  real mean_serial_interval_real;

  // reverse infection and test delay profiles and accumulated result delay profile
  {
    real mean_serial_interval = 0.0;
    for (i in 1:Tip) {
      infprofile_rev[Tip-i+1] = infprofile[i];
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
    for (i in 1:Ttdp)
      testdelayprofile_rev[i] = testdelayprofile[Ttdp-i+1];
  }
}

parameters {
  real mu; // prior mean for log(Rt)
  real<lower=0> sigma; // prior scale for log(Rt)
  real<lower=0,upper=1> alpha1; // 1-autocorrelation for log(Rt)
  real<lower=0> phi_latent; // dispersion for negative binomial latent process
  real<lower=0> phi_observed; // dispersion for negative binomial observation
  vector[Nstep+Nproj] Reta;
  vector[Tlik] Ceta;
  real<lower=0> xi;
  simplex[Trdp] resultdelayprofile;
}

transformed parameters {
  real alpha;
  vector[Nstep+Nproj] Rt;
  vector[Tcur] Xt;
  vector[Tlik] Ecount;
  vector[Trdp] resultdelayprofile_revcum;


  if (fixed_gp_time_length_scale <= 0.0) {
    alpha = pow(1.0-alpha1, Tstep/gp_time_scale); // autocorrelation for log(Rt)
  } else {
    alpha = exp(-Tstep/fixed_gp_time_length_scale);
  }
  { // AR1 process for log(Rt)
    vector[Nstep+Nproj] Zt = rep_vector(0.0,Nstep+Nproj);
    real delta = sqrt(1.0-alpha);
    Zt[1] += Reta[1];
    for (i in 2:Nstep+Nproj)
      Zt[i] += alpha * Zt[i-1] + delta * Reta[i];
    Rt = exp(mu + sigma * Zt);
  }
  { // results delay distribution
    real s = 0.0;
    for (i in 1:Trdp) {
      s += resultdelayprofile[i];
      resultdelayprofile_revcum[Trdp-i+1] = s;
    }
  }
  { // latent renewal process 
    // negative binomial noise approximated by gaussian
    for (t in 1:Tcond) {
      Xt[t] = xi + 
        (1.0-mean_serial_interval_real) * Count[1,t+mean_serial_interval_int] + 
        mean_serial_interval_real * Count[1,t+mean_serial_interval_int+1];
    }
    for (i in 1:Nstep) {
      for (j in 1:Tstep) {
        int s = (i-1)*Tstep + j;
        int t = Tcond + s;
        int L = min(Tip,t-1);

        real Einfection = Rt[i] * (xi + dot_product(
            Xt[t-L:t-1], 
            infprofile_rev[Tip-L+1:Tip]
        ));
        //approximate poisson with log normal with same mean/variance
        Xt[t] = fabs(
          Einfection + 
          //sqrt(Einfection) * Ceta[s] // Poisson
          sqrt((1.0+phi_latent) * Einfection) * Ceta[s] // neg-binomial
        );

        Ecount[s] = dot_product(
            Xt[t-Ttdp+1:t], 
            testdelayprofile_rev
        );
    } } 
  }
}

model {
  mu ~ normal(0.0, mu_scale);
  sigma ~ normal(0.0, sigma_scale);
  alpha1 ~ normal(0.0, gp_time_decay_scale);
  phi_latent ~ normal(0.0, phi_latent_scale);
  phi_observed ~ normal(0.0, phi_observed_scale);
  xi ~ normal(0.0, xi_scale);
  Reta ~ std_normal();
  Ceta ~ std_normal();

  {
    resultdelayprofile ~ dirichlet(resultdelayalpha[1:Trdp]);
  }

  {
    for (t in Tcond+1:Tcur-Trdp+1) {
      int s = t-Tcond;
      Count[1,t] ~ neg_binomial_2(Ecount[s], 1.0 / phi_observed);
    }
    for (i in 2:Trdp) {
      int t = Tcur-Trdp+i;
      int s = t-Tcond;
      real ec = Ecount[s] * resultdelayprofile_revcum[i];
      Count[1,t] ~ neg_binomial_2(ec, 1.0 / phi_observed);
    }
  }
}

generated quantities {
  int Crecon[Tcur];
  vector[Tstep*Nstep] Cpred;
  int Noutliers = 0;
  real meandelay = 0.0;
  real denomdelay = 0.0;
  vector[Tstep*Nproj] Xt_proj;
  vector[Tstep*Nproj] Count_proj;
  vector[Tpred] Ppred;
  real gp_time_length_scale = -Tstep/log(alpha);

  for (t in 1:(Tstep*Nstep)) {
    Cpred[t] = 0.0;
  }

  for (t in 1:Tcond)
    Crecon[t] = Count[1,t];
  for (t in Tcond+1:Tcond+Tlik) 
    Crecon[t] = 0;
  for (t in Tcond+1:Tcur) {
    int s = t-Tcond;
    int c = Count[1,t];
    real Ec = Ecount[s];
    real psi = Ec / phi_observed;
    vector[Ttdp] Precon;
    int Crecon_t[Ttdp];
    if (c>outlier_count_threshold && 
               neg_binomial_2_cdf(c,Ec,psi)>outlier_prob_threshold) {
      c = max(outlier_count_threshold,neg_binomial_2_rng(Ec,psi));
      Noutliers += 1;
    }
    Precon = testdelayprofile_rev .* Xt[t-Ttdp+1:t];
    Cpred[s] = sum(Precon);
    if (c==0) continue;
    denomdelay += c;
    if (reconstruct_infections) {
      Precon /= sum(Precon);
      Crecon_t = multinomial_rng(Precon,c);
      for (i in 1:Ttdp) {
        Crecon[t-Ttdp+i] += Crecon_t[i];
        meandelay += c*Precon[i]*(Ttdp-i);
      }
    } else {
      Crecon[t] = c;
    }
  }
  meandelay /= denomdelay;

  {
    vector[Tcur+Tstep*Nproj] xlatent;
    xlatent[1:Tcur] = Xt;
    for (i in Nstep+1:Nstep+Nproj) {
      for (j in 1:Tstep) {
        int s = (i-1)*Tstep + j;
        int t = Tcond + s;
        int L = min(Tip,t-1);
        real Einfection = Rt[i] * (xi + dot_product(
            xlatent[t-L:t-1], 
            infprofile_rev[Tip-L+1:Tip]
        ));
        //approximate poisson with log normal with same mean/variance
        xlatent[t] = fabs(
          Einfection + 
          //sqrt(Einfection) * normal_rng(0.0,1.0) // Poisson
          sqrt((1.0+phi_latent) * Einfection) * normal_rng(0.0,1.0) // NB
        );

        Count_proj[t-Tcur] = dot_product(
            xlatent[t-Ttdp+1:t], 
            testdelayprofile_rev
        );
    } } 
    Xt_proj = xlatent[Tcur+1:Tcur+Tstep*Nproj];

    for (m in 1:Nproj) {
      for (t in 1:Tstep) {
        int i = (m-1) * Tstep + t;
        if (i > Tpred) {
          break;
        }
        Ppred[i] = exp(neg_binomial_2_lpmf(Count[1, Tcur + i] |
            Count_proj[i],
            1.0 / phi_observed
        ));
      }
    }
  }


} 

