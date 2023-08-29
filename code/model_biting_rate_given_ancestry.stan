// Probability of biting a human given a bite

data {
 int omega_ancestry_N;                          // number of observations
 vector[omega_ancestry_N] omega_ancestry;       // vector of prob. biting human | ancestry
 vector[omega_ancestry_N] aa_ancestry;          // vector of proportion aa ancestry
 
 // new ancestry data for predictions
 int ancestry_N_new;                            // number of new ancestry values
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
  omega_ancestry_d ~ uniform(0,1);             // prior for d
  omega_ancestry_e ~ uniform(0,10);             // prior for e
  omega_ancestry_sigma ~ uniform(0,100);        // prior for sigma
  
  // omega prob biting human | ancestry
  for(h in 1:omega_ancestry_N){
    real omega_ancestry_mu = omega_ancestry_constant + ((omega_ancestry_d - omega_ancestry_constant)/(1 + (omega_ancestry_e / aa_ancestry[h])));
    omega_ancestry[h] ~ normal(omega_ancestry_mu, omega_ancestry_sigma);
  }
} 

generated quantities {
  // posterior predictive check values  
  real omega_ancestry_ppc[omega_ancestry_N];
  vector[ancestry_N_new] omega_ancestry_new;
  
  // new vectors for parameters
  vector[ancestry_N_new] alpha;                 // biting rate
  vector[ancestry_N_new] b;                     // probability of mosquito infectiousness
  vector[ancestry_N_new] lf;                    // lifepsan
  vector[ancestry_N_new] EIR;                   // extrinsic incubation rate
  vector[ancestry_N_new] pMI;                   // probability of mosquito infection
  vector[ancestry_N_new] NmNh;                  // ratio of mosquitoes to humans
  vector[ancestry_N_new] delta;                 // intrinsic incubation period
  vector[ancestry_N_new] mu_h;                  // human mortality rate
  vector[ancestry_N_new] gamma;                 // Human infectivity period
  vector[ancestry_N_new] R0_ancestry;           // R0

  // omega ppc prob biting human | ancestry
  for (o in 1:omega_ancestry_N){
    real omega_ancestry_mu_ppc = omega_ancestry_constant + ((omega_ancestry_d - omega_ancestry_constant)/(1 + (omega_ancestry_e / aa_ancestry[o])));
    omega_ancestry_ppc[o] = normal_rng(omega_ancestry_mu_ppc, omega_ancestry_sigma);
  }
  
  // omega prob biting human | ancestry across new gradient
  for (yy in 1:ancestry_N_new){
    omega_ancestry_new[yy] = normal_rng(omega_ancestry_constant + ((omega_ancestry_d - omega_ancestry_constant)/(1 + (omega_ancestry_e / aa_ancestry_new[yy]))), omega_ancestry_sigma);

    alpha[yy] = normal_rng(2.17, 0.62);
    b[yy] = normal_rng(0.363, 0.217);
    EIR[yy] = normal_rng(8.97, 1.22);
    lf[yy] = normal_rng(16.5, 2.4);
    pMI[yy] = normal_rng(0.67, 0.17);
    NmNh[yy] = normal_rng(2,0.5);
    mu_h[yy] = normal_rng(4.25E-05, 0.00005);
    delta[yy] = 1/gamma_rng(5.9, 0.5);                  //zikv from https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0004726
    gamma[yy] = 1/gamma_rng(5.0, 0.5);

   // R0
    R0_ancestry[yy] = sqrt(
      (omega_ancestry_new[yy] * alpha[yy] * b[yy] * (EIR[yy] / ((1/lf[yy]) * ((1/lf[yy]) + EIR[yy])))) *
      (alpha[yy] * pMI[yy] * NmNh[yy] * (delta[yy] / ((delta[yy] + mu_h[yy]) * (gamma[yy] + mu_h[yy]))))
      );

  }
  
}
