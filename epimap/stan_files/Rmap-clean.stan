functions {
  real adjustedR(real RZ, real dispersion) {
    return RZ 
      - 2.0 * RZ * normal_cdf(-sqrt(RZ/(1.0+dispersion)),0.0,1.0) 
      + sqrt(2.0*(1.0+dispersion)*RZ/pi())*exp(-RZ/(2.0*(1.0+dispersion)));
  }
}

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
  real mu; // prior mean for log(Rx)
  real<lower=0> sigma; // prior scale for log(Rx)
  real<lower=0,upper=1> alpha1; // 1-autocorrelation for log(Rx)
  real<lower=0> phi_latent; // dispersion for negative binomial latent process
  real<lower=0> phi_observed; // dispersion for negative binomial observation
  vector[Nstep+Nproj] Reta;
  vector[Tlik] Ceta;
  real<lower=0> xi;
  simplex[Trdp] resultdelayprofile;
  simplex[7] weekly_case_variations;
}

transformed parameters {
  real alpha;
  vector[Nstep+Nproj] Rx;
  vector[Tcur] Xt;
  vector[Tlik] Ecount;
  vector[Trdp] resultdelayprofile_revcum;


  if (fixed_gp_time_length_scale <= 0.0) {
    alpha = pow(1.0-alpha1, Tstep/gp_time_scale); // autocorrelation for log(Rx)
  } else {
    alpha = exp(-Tstep/fixed_gp_time_length_scale);
  }
  { // AR1 process for log(Rx)
    vector[Nstep+Nproj] Zt = rep_vector(0.0,Nstep+Nproj);
    real delta = sqrt(1.0-alpha);
    Zt[1] += Reta[1];
    for (i in 2:Nstep+Nproj)
      Zt[i] += alpha * Zt[i-1] + delta * Reta[i];
    //Rx = exp(mu + sigma * Zt); 
    Rx = exp(sigma * Zt); // don't want drift when forecasting
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

        real Einfection = Rx[i] * (xi + dot_product(
            Xt[t-L:t-1], 
            infprofile_rev[Tip-L+1:Tip]
        ));
        //approximate poisson with log normal with same mean/variance
        Xt[t] = fabs(
          Einfection + 
          sqrt(Einfection) * Ceta[s] // Poisson
          //sqrt((1.0+1.0/phi_latent) * Einfection) * Ceta[s] // neg-binomial
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
  weekly_case_variations ~ dirichlet(rep_vector(7.0,7));

  {
    resultdelayprofile ~ dirichlet(resultdelayalpha[1:Trdp]);
  }

  {
    for (t in Tcond+1:Tcur-Trdp+1) {
      int s = t-Tcond;
      int d = (t % 7)+1;
      real ec = (7.0*weekly_case_variations[d]) * Ecount[s];
      Count[1,t] ~ neg_binomial_2(fmax(1e-3, ec), fmax(1e-3,ec / phi_observed));
    }
    for (i in 2:Trdp) {
      int t = Tcur-Trdp+i;
      int s = t-Tcond;
      int d = (t % 7)+1;
      real ec = (7.0*weekly_case_variations[d]) *
        Ecount[s] * resultdelayprofile_revcum[i];
      Count[1,t] ~ neg_binomial_2(
        fmax(1e-3, ec), 
        fmax(1e-3, ec / phi_observed)
      );
    }
  }
}

generated quantities {
  vector[Nstep+Nproj] Rt;
  int Crecon[Tcur];
  int Noutliers = 0;
  real meandelay = 0.0;
  real denomdelay = 0.0;
  vector[Tcur+Tstep*Nproj] Xt_proj;
  vector[Tstep*Nstep] Cpred;
  vector[Tstep*Nproj] Cproj;
  vector[Tstep*Nstep] Xpred;
  vector[Tstep*Nproj] Xproj;
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
    real ecpred;
    real psi = Ec / phi_observed;
    vector[Ttdp] Precon;
    int Crecon_t[Ttdp];
    int d = (t % 7)+1;
    if (c>outlier_count_threshold && 
               neg_binomial_2_cdf(c,Ec,psi)>outlier_prob_threshold) {
      c = max(outlier_count_threshold,neg_binomial_2_rng(Ec,psi));
      Noutliers += 1;
    }
    Precon = testdelayprofile_rev .* Xt[t-Ttdp+1:t];
    ecpred =  (7.0*weekly_case_variations[d]) * sum(Precon);
    //if (ecpred<1e-3 || ecpred>1e3 || phi_observed<1e-3 || phi_observed>1e3) {
    //  print("ecpred ",ecpred, " phi_observed ",phi_observed);
    //}
    Cpred[s] = neg_binomial_2_rng(
      fmin(1e5, fmax(1e-3, ecpred)), 
      fmax(1e-3, ecpred / phi_observed)
    );

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
        int d = (t % 7)+1;
        real ecproj;
        real Einfection = Rx[i] * (xi + dot_product(
            xlatent[t-L:t-1], 
            infprofile_rev[Tip-L+1:Tip]
        ));
        //approximate poisson with log normal with same mean/variance
        xlatent[t] = fabs(
          Einfection + 
          sqrt(Einfection) * normal_rng(0.0,1.0) // Poisson
          //sqrt((1.0+1.0/phi_latent) * Einfection) * normal_rng(0.0,1.0) // NB
        );

        ecproj = (7.0*weekly_case_variations[d]) * 
          dot_product(xlatent[t-Ttdp+1:t], testdelayprofile_rev);
        Cproj[t-Tcur] = neg_binomial_2_rng(
          fmin(1e5, fmax(1e-3, ecproj)), 
          fmax(1e-3, ecproj / phi_observed)
        );
    } } 
    Xt_proj = xlatent;
    Xpred = xlatent[Tcond+1:Tcur];
    Xproj = xlatent[Tcur+1:Tcur+Tstep*Nproj];

    for (m in 1:Nproj) {
      for (t in 1:Tstep) {
        int i = (m-1) * Tstep + t;
        if (i > Tpred) {
          break;
        }
        Ppred[i] = exp(neg_binomial_2_lpmf(Count[1, Tcur + i] |
            fmax(1e-3,Cproj[i]),
            fmax(1e-3,Cproj[i]) / phi_observed
        ));
      }
    }

    for (i in 1:Nstep+Nproj) {
      real sumZ = 0.0;
      real sumRZ = 0.0;
      for (s in 1:Tstep) {
        int t = Tcond + (i-1)*Tstep + s;
        int L = min(Tip,t-1);
        real Z = (xi + dot_product(
            Xt_proj[t-L:t-1], 
            infprofile_rev[Tip-L+1:Tip]
        ));
        real RZ = Rx[i] * Z;
        sumZ += Z;
        sumRZ += RZ;
        // sumRZ += Xt_proj[t];
        // sumRZ += adjustedR(RZ,1.0/phi_latent);
      }
      Rt[i] = sumRZ / sumZ;
    }

  }

} 

