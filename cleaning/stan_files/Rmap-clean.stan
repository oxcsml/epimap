data {
  int<lower=1> Tall;  //
  int<lower=1> Tstep; // length of step (7 days = 1 week)
  int<lower=1> Nstep; // number of steps
  int Count[1,Tall];
  int<lower=1> Tip;
  vector[Tip] infprofile;
  int<lower=1> Ttdp;
  vector[Ttdp] testdelayprofile;
  int<lower=1> Trdp;
  real<lower=0,upper=1> resultdelaydecay;
  real<lower=1> resultdelaystrength;
  real<lower=0> mu_scale;
  real<lower=0> sigma_scale;
  real<lower=0> alpha_scale;
  real<lower=0> phi_latent_scale;
  real<lower=0> phi_observed_scale;
  real<lower=0,upper=1> outlier_prob_threshold;
  int<lower=1> outlier_count_threshold;
  real<lower=0> xi_scale;
  int<lower=0,upper=1> reconstruct_infections;
}

transformed data {
  int Tlik = Nstep*Tstep;
  int Tcond = Tall-Tlik;
  vector[Tip] infprofile_rev;
  vector[Ttdp] testdelayprofile_rev;

  // reverse infection and test delay profiles and accumulated result delay profile
  {
    for (i in 1:Tip)
      infprofile_rev[i] = infprofile[Tip-i+1];
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
  vector[Nstep] Reta;
  vector[Tlik] Ceta;
  real<lower=0> xi;
  simplex[Trdp] resultdelayprofile;
}

transformed parameters {
  real alpha = 1-alpha1; // autocorrelation for log(Rt)
  vector[Nstep] Rt;
  vector[Tall] Clatent;
  vector[Tlik] Ecount;
  vector[Trdp] resultdelayprofile_revcum;

  { // AR1 process for log(Rt)
    vector[Nstep] Zt = rep_vector(0.0,Nstep);
    real delta = sqrt(1.0-alpha);
    Zt[1] += Reta[1];
    for (i in 2:Nstep)
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
    for (t in 1:Tcond) 
      Clatent[t] = xi + Count[1,t];
    for (i in 1:Nstep) {
      for (j in 1:Tstep) {
        int s = (i-1)*Tstep + j;
        int t = Tcond + s;
        int L = min(Tip,t-1);
        real Einfection = Rt[i] * (xi + dot_product(
            Clatent[t-L:t-1], 
            infprofile_rev[Tip-L+1:Tip]
        ));
        Clatent[t] = fabs(
          Einfection + 
          sqrt(Einfection) * Ceta[s]
          // sqrt((1.0+phi_latent) * Einfection) * Ceta[s]
        );

        //approximate poisson with log normal with same mean/variance
        //real logEinfection = log(Einfection);
        //real logEinfection1 = log(Einfection+1.0);
        //Clatent[t] = exp(
        //    1.5 * logEinfection - .5 * logEinfection1 +
        //    sqrt(logEinfection1 - logEinfection) * Ceta[s]
        //);

        Ecount[s] = dot_product(
            Clatent[t-Ttdp+1:t], 
            testdelayprofile_rev
        );
    } } 
  }
}

model {
  mu ~ normal(0.0, mu_scale);
  sigma ~ normal(0.0, sigma_scale);
  alpha1 ~ normal(0.0, alpha_scale);
  phi_latent ~ normal(0.0, phi_latent_scale);
  phi_observed ~ normal(0.0, phi_observed_scale);
  xi ~ normal(0.0, xi_scale);
  Reta ~ std_normal();
  Ceta ~ std_normal();

  {
    vector[Trdp] diralpha;
    for (i in 1:Trdp) 
      diralpha[i] = resultdelaystrength * resultdelaydecay^i;
    resultdelayprofile ~ dirichlet(diralpha);
  }

  {
    for (t in Tcond+1:Tall-Trdp+1) {
      int s = t-Tcond;
      Count[1,t] ~ neg_binomial_2(Ecount[s], Ecount[s] / phi_observed);
    }
    for (i in 2:Trdp) {
      int t = Tall-Trdp+i;
      int s = t-Tcond;
      real ec = Ecount[s] * resultdelayprofile_revcum[i];
      Count[1,t] ~ neg_binomial_2(ec, ec / phi_observed);
    }
  }
}

generated quantities {
  int Crecon[Tall];
  int Noutliers = 0;
  real meandelay = 0.0;
  real denomdelay = 0.0;

  for (t in 1:Tcond)
    Crecon[t] = Count[1,t];
  for (t in Tcond+1:Tcond+Tlik) 
    Crecon[t] = 0;
  for (t in Tcond+1:Tall) {
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
    if (c==0) continue;
    denomdelay += c;
    if (reconstruct_infections) {
      Precon = testdelayprofile_rev .* Clatent[t-Ttdp+1:t];
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
} 

