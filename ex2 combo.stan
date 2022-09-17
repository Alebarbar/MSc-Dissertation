data {
  int<lower=1> I;
  int<lower=1> J;
  real<lower=0> sigma_a_upper;
  real<lower=0> sigma_b_upper;
  real<lower=0,upper=1> eta_a;
  real<lower=0,upper=1> eta_b;
  real<lower=0,upper=1> zeta;
  int<lower=0,upper=1> rho_prior;
  int<lower=0,upper=1> bias_prior;
  matrix<lower=0,upper=1>[I,J] y;
}

transformed data {
  matrix<lower=0,upper=1>[I,J] x;
  for (j in 1:J){
    for (i in 1:I){
      if (y[i,j] > 1-10^(-16)){
        x[i,j] = 1-10^(-16);
      }else if (y[i,j] <10^(-16)){
        x[i,j] = 10^(-16);
      }else{
        x[i,j] = y[i,j];
       }
    }
  } 
}

parameters {
  real<lower=0,upper=sigma_a_upper> sigma_a;
  real<lower=0,upper=sigma_b_upper> sigma_b;
  matrix[J,2] effects;
  real<lower=0,upper=1> rho;
  real<lower=0,upper=1> bias;
}

transformed parameters {
  vector<lower=0>[J] alpha;
  if(rho_prior==1){
    alpha = exp(effects[,1]*sigma_a);
  }else{
    alpha = exp(effects[,1]);
  }
  vector<lower=0>[J] beta;
  if(rho_prior==1){
    beta = exp(effects[,1]*rho*sigma_b + effects[,2]*sqrt(1-rho^2)*sigma_b);
  }else{
    beta = exp(effects[,2]);
  }

}

model {
  sigma_a ~ uniform(0,sigma_a_upper);
  sigma_b ~ uniform(0,sigma_b_upper);
  rho ~ uniform(0,1);
  bias ~ uniform(0,1);
  if(rho_prior==1){
    for (k in 1:2){
      for (j in 1:J){
        effects[j,k] ~ normal(0,1);
      }
    }
  }else{
    target += eta_a*normal_lpdf(effects[,1]|rep_vector(0,J),sigma_a);
    target += eta_b*normal_lpdf(effects[,2]|rep_vector(0,J),sigma_b);
  }
  
  for (j in 1:J){
    if (bias_prior==1){
      for (i in 1:I){
        target += log_sum_exp(log(1-bias)+beta_lpdf(x[i,j]|alpha[j],beta[j]),log(bias));
      }
    }else{
      target += zeta*beta_lpdf(x[,j]|alpha[j],beta[j]);
    }
  }
}

