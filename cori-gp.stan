// saved as schools.stan
data {
  int<lower=1> N;                 // number of regions
  int<lower=1> T;                 // number of days
  int<lower=1> T0;                // number of recent days for likelihood computation
  int<lower=1> D;                 // length of infection profile. 
  int C[N, T];                    // case counts
  vector[D] infprofile;           // infection profile
}

transformed data {
  vector[1] x[N];                 // !!TEMP!! region covariate
  vector[N] mu;                   // GP mean 
  vector[D] infprofile_rev;       // reversed infection profile
  vector[T] Cvec[N];              // counts in real vector form
  vector[T0] convout[N];          // convolution between C and infprofile

  for (i in 1:N) 
    x[i][1] = i;
  mu = rep_vector(0,N); 

  for (i in 1:D)
    infprofile_rev[i] = infprofile[D-i+1];
  for (j in 1:N) {
    for (i in 1:T)
      Cvec[j][i] = C[j,i];
    for (i in 1:T0) {
      int L = min(D,T-T0+i-1); // length of infection profile that overlaps with case counts in past days
      convout[j][i] = dot_product(Cvec[j][T-T0-L+i:T-T0+i-1], infprofile_rev[D-L+1:D]);
    }
  }
}

parameters {
  real<lower=0> length_scale;
  real<lower=0> func_sigma;
  real<lower=0> data_sigma;
  vector[N] eta;
}

transformed parameters {
  vector[N] Rt;                 // instantaneous reproduction number

  {
    real data_sigma2;
    matrix[N,N] K;
    matrix[N,N] L;

    data_sigma2 = square(data_sigma);

    K = cov_exp_quad(x, func_sigma, length_scale); // kernel
    for (i in 1:N) 
      K[i,i] = K[i,i] + data_sigma2;
  
    L = cholesky_decompose(K);
    Rt = exp(L * eta); 
  }
}
 
model {

  // GP
  length_scale ~ inv_gamma(5,5);
  func_sigma ~ std_normal();
  data_sigma ~ std_normal();
  eta ~ std_normal();

  for (j in 1:N) {
    for (i in 1:T0) {
      C[j][T-T0+i] ~ poisson(Rt[j] * convout[j][i]);
    }
  }
}


