
data {
  int<lower=1> I;
  int<lower=1> J;
  int<lower=1> K;
  int<lower=1> D;
  real<lower=0> sigma_v_upper;
  real<lower=0> sigma_u_upper;
  real<lower=0> beta_sd;
  real<lower=0,upper=1> rho;
  real<lower=0,upper=1> eta_v;
  real<lower=0,upper=1> eta_u;
  real<lower=0,upper=1> zeta;
  array[I,J,K] int<lower=0,upper=1> y;
}


parameters {
  matrix[J+1,K] effects;
  array[I,J,K] real Z;
  real<lower=0> sigma_v;
  real<lower=0> sigma_u;
  real beta_0;
  real<lower=0> beta_sq;
}

transformed parameters {
  vector[K] v;
  matrix[J,K] u;
  for (k in 1:K){
    v[k] = sigma_v*((sqrt(1+J*rho)-sqrt(1-rho))*sum(effects[,k])/(J+1)+sqrt(1-rho)*effects[1,k]);
    for (j in 1:J){
      u[j,k]= sigma_u*((sqrt(1+J*rho)-sqrt(1-rho))*sum(effects[,k])/(J+1)+sqrt(1-rho)*effects[j+1,k]);
    }
  }
}

model {
  sigma_v ~ uniform(0,sigma_v_upper);
  sigma_u ~ uniform(0,sigma_u_upper);
  for (k in 1:K){
    target += eta_v*normal_lpdf(effects[1,k]|0,1);
    for (j in 1:J){
      target += eta_u*normal_lpdf(effects[j+1,k]|0,1);
    }
  }
  beta_0 ~ normal(0,beta_sd);
  (beta_sq/beta_sd)^2 ~ chi_square(D);
  array[I,J,K] real p;
  for (k in 1:K){
    for (j in 1:J){
      for (i in 1:I){
        Z[i,j,k] ~ normal(0,1);
        p[i,j,k] = beta_0 + Z[i,j,k] * beta_sq + u[j,k] + v[k];
        target += zeta*bernoulli_logit_lpmf(y[i,j,k]|p[i,j,k]);
      }
    }
  }
}

generated quantities{
  vector[D] beta;
  beta = rep_vector(beta_sq/sqrt(D),D);
}
