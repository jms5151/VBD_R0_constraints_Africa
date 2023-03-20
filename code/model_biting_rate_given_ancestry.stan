// maybe want to rename, this isn't really alpha (biting rate)
// it's probability of biting a human given a bite

data {
 int alpha_ancestry_N;                          // number of observations
 vector[alpha_ancestry_N] alpha_ancestry;       // vector of prob. biting | ancestry
 vector[alpha_ancestry_N] aa_ancestry;          // vector of proportion aa ancestry
 
 int ancestry_N_new;                            // number of new observations to predict on
 vector[ancestry_N_new] aa_ancestry_new;        // vector of proportion aa ancestry 
}

parameters {
 real<lower=0> alpha_ancestry_constant;         // parameter c, lower limit
 real<lower=0> alpha_ancestry_d;                // parameter d, upper limit
 real<lower=0> alpha_ancestry_e;                // parameter e, dose responding to halfway between c and d
 real<lower=0> alpha_ancestry_sigma;            // noise
}

model {
  
  alpha_ancestry_constant ~ uniform(0,1);       // prior for c
  alpha_ancestry_d ~ uniform(0,10);             // prior for d
  alpha_ancestry_e ~ uniform(0,10);             // prior for e
  alpha_ancestry_sigma ~ uniform(0,100);        // prior for sigma
  
  // alpha (maybe something new?) prob biting | ancestry
  for(h in 1:alpha_ancestry_N){
    real alpha_ancestry_mu = alpha_ancestry_constant + ((alpha_ancestry_d - alpha_ancestry_constant)/(1 + (alpha_ancestry_e / aa_ancestry[h])));
    alpha_ancestry[h] ~ normal(alpha_ancestry_mu, alpha_ancestry_sigma);
  }
} 

generated quantities {
  
  real alpha_ancestry_ppc[alpha_ancestry_N];
  vector[ancestry_N_new] alpha_ancestry_new;
  
  // alpha ppc (maybe something new?) prob biting | ancestry
  for (o in 1:alpha_ancestry_N){
    real alpha_ancestry_mu_ppc = alpha_ancestry_constant + ((alpha_ancestry_d - alpha_ancestry_constant)/(1 + (alpha_ancestry_e / aa_ancestry[o])));
    alpha_ancestry_ppc[o] = normal_rng(alpha_ancestry_mu_ppc, alpha_ancestry_sigma);
  }
  
  // alpha (maybe something new?) prob biting | ancestry across new gradient
  for (yy in 1:ancestry_N_new){
    alpha_ancestry_new[yy] = normal_rng(alpha_ancestry_constant + ((alpha_ancestry_d - alpha_ancestry_constant)/(1 + (alpha_ancestry_e / aa_ancestry_new[yy]))), alpha_ancestry_sigma);
  }
}