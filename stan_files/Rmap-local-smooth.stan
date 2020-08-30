data {
  int<lower=1> Tall;  //
  int<lower=1> Tstep; // length of step (7 days = 1 week)
  int<lower=1> Nstep; // number of steps
  int Count[1,Tall];
  int<lower=1> Tip;
  vector[Tip] infprofile;
  int<lower=1> Tdp;
  vector[Tdp] delayprofile;
  real<lower=0> mu_scale;
  real<lower=0> sigma_scale;
  real<lower=0> alpha_scale;
  real<lower=0> phi_latent_scale;
  real<lower=0> phi_observed_scale;
}

transformed data {
  int Tlik = Nstep*Tstep;
  int Tcond = Tall-Tlik;
  vector[Tip] infprofile_rev;
  vector[Tdp] delayprofile_rev;

  // reverse infection and delay profile
  for (i in 1:Tip)
    infprofile_rev[i] = infprofile[Tip-i+1];

  for (i in 1:Tdp)
    delayprofile_rev[i] = delayprofile[Tdp-i+1];
}

parameters {
  real mu; // prior mean for log(Rt)
  real<lower=0> sigma; // prior scale for log(Rt)
  real<lower=0,upper=1> alpha1; // 1-autocorrelation for log(Rt)
  real<lower=0> phi_latent; // dispersion for negative binomial latent process
  real<lower=0> phi_observed; // dispersion for negative binomial observation
  vector[Nstep] Reta;
  vector[Tlik] Ceta;
}

transformed parameters {
  real alpha = 1-alpha1; // autocorrelation for log(Rt)
  vector[Nstep] Rt;
  vector[Tlik] Cinfer;

  { // AR1 process for log(Rt)
    vector[Nstep] Zt;
    real delta = sqrt(1.0-alpha);
    Zt[1] = sigma * Reta[1];
    for (i in 2:Nstep) {
      Zt[i] = alpha * Zt[i-1] + delta * sigma * Reta[i];
    }
    Rt = exp(mu + Zt);
  }
  { // latent renewal process 
    // negative binomial noise approximated by gaussian
    vector[Tall] Clatent;
    for (t in 1:Tcond) 
      Clatent[t] = Count[1,t];
    for (i in 1:Nstep) {
      for (j in 1:Tstep) {
        int s = (i-1)*Tstep + j;
        int t = Tcond + s;
        int L = min(Tip,t-1);
        real Einf = Rt[i] * dot_product(
            Clatent[t-L:t-1], 
            infprofile_rev[Tip-L+1:Tip]
        );
        Clatent[t] = fabs(Einf + sqrt((1.0+phi_latent) * Einf) * Ceta[s]);
    } } 
    Cinfer = Clatent[Tcond+1:Tall];
  }
}

model {
  mu ~ normal(0.0, mu_scale);
  sigma ~ normal(0.0, sigma_scale);
  alpha1 ~ normal(0.0, alpha_scale);
  phi_latent ~ normal(0.0, phi_latent_scale);
  phi_observed ~ normal(0.0, phi_observed_scale);
  Reta ~ std_normal();
  Ceta ~ std_normal();

  {
    vector[Tall] Clatent;
    for (t in 1:Tcond) 
      Clatent[t] = Count[1,t];
    Clatent[Tcond+1:Tall] = Cinfer;
    for (t in Tcond+1:Tall) {
      real Ecase = dot_product(Clatent[t-Tdp:t-1], delayprofile_rev);
      Count[1,t] ~ neg_binomial_2(Ecase, Ecase / phi_observed);
    }
  }
}

generated quantities {
} 

