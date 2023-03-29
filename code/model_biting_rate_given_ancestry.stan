// it's probability of biting a human given a bite

data {
 int omega_ancestry_N;                          // number of observations
 vector[omega_ancestry_N] omega_ancestry;       // vector of prob. biting human | ancestry
 vector[omega_ancestry_N] aa_ancestry;          // vector of proportion aa ancestry
 
 int ancestry_N_new;                            // number of new observations to predict on
 vector[ancestry_N_new] aa_ancestry_new;        // vector of proportion aa ancestry 
}

parameters {
 real<lower=0> omega_ancestry_constant;         // parameter c, lower limit
 real<lower=0, upper=1> omega_ancestry_d;       // parameter d, upper limit
 real<lower=0> omega_ancestry_e;                // parameter e, dose responding to halfway between c and d
 real<lower=0> omega_ancestry_sigma;            // noise
}

model {
  
  omega_ancestry_constant ~ uniform(0,1);       // prior for c
  omega_ancestry_d ~ uniform(0,10);             // prior for d
  omega_ancestry_e ~ uniform(0,10);             // prior for e
  omega_ancestry_sigma ~ uniform(0,100);        // prior for sigma
  
  // omega prob biting human | ancestry
  for(h in 1:omega_ancestry_N){
    real omega_ancestry_mu = omega_ancestry_constant + ((omega_ancestry_d - omega_ancestry_constant)/(1 + (omega_ancestry_e / aa_ancestry[h])));
    omega_ancestry[h] ~ normal(omega_ancestry_mu, omega_ancestry_sigma);
  }
} 

generated quantities {
  
  real omega_ancestry_ppc[omega_ancestry_N];
  vector[ancestry_N_new] omega_ancestry_new;
  
  // omega ppc prob biting human | ancestry
  for (o in 1:omega_ancestry_N){
    real omega_ancestry_mu_ppc = omega_ancestry_constant + ((omega_ancestry_d - omega_ancestry_constant)/(1 + (omega_ancestry_e / aa_ancestry[o])));
    omega_ancestry_ppc[o] = normal_rng(omega_ancestry_mu_ppc, omega_ancestry_sigma);
  }
  
  // omega prob biting human | ancestry across new gradient
  for (yy in 1:ancestry_N_new){
    omega_ancestry_new[yy] = normal_rng(omega_ancestry_constant + ((omega_ancestry_d - omega_ancestry_constant)/(1 + (omega_ancestry_e / aa_ancestry_new[yy]))), omega_ancestry_sigma);
  }
}