data {
  
  // omega (probability of biting a human given a bite)
  int omega_ancestry_N;                                         // number of observations
  vector[omega_ancestry_N] omega_ancestry;                      // vector of prob. biting human | ancestry
  vector[omega_ancestry_N] omega_ancestry_aa;                   // vector of proportion aa ancestry

  // alpha (biting rate)
  int alpha_climate_N;                                          // number of observations
  vector[alpha_climate_N] alpha_climate_temp;                   // vector of temperatures
  vector[alpha_climate_N] alpha_climate;                        // vector of trait alpha

  // b (prob mosquito infectiousness)
  int b_climate_N;                                              // number of observations
  vector[b_climate_N] b_climate_temp;                           // vector of temperatures
  vector[b_climate_N] b_climate;                                // vector of trait b
  
  // pMI climate (prob mosquito infection)
  int pMI_climate_N;                                            // number of observations
  vector[pMI_climate_N] pMI_climate_temp;                       // vector of temperatures
  vector[pMI_climate_N] pMI_climate;                            // vector of trait pMI
  
  // pMI ancestry (prob mosquito infection)
  int<lower=0> pMI_ancestry_N;                                  // number of observations
  int<lower=0> pMI_ancestry_Trials[pMI_ancestry_N];             // number of trials for each observation
  int<lower=0> pMI_ancestry_Infected[pMI_ancestry_N];           // number of successes for each observation
  matrix[pMI_ancestry_N, 2] pMI_ancestry_X;                     // matrix of covariates (continuous)
  int<lower=0, upper=1> pMI_ancestry_Virus[pMI_ancestry_N];     // categorical variable

  // EIR (extrinsic incubation rate)
  int EIR_climate_N;                                            // number of observations
  vector[EIR_climate_N] EIR_climate_temp;                       // vector of temperatures
  vector[EIR_climate_N] EIR_climate;                            // vector of trait EIR

  // lifespan (1/mosquito mortality rate)
  int lf_climate_N;                                             // number of observations
  vector[lf_climate_N] lf_climate_temp;                         // vector of temperatures
  vector[lf_climate_N] lf_climate;                              // vector of trait lf
  
  // new data
  int N_new;
  vector[N_new] temp_new;
  vector[N_new] aa_new;
  int<lower=0, upper=1> Virus_new[N_new];
  matrix[N_new, 2] X_new;

  // // new climate data for predictions 
  // int climate_N_new;                                            // number of new temperatures
  // vector[climate_N_new] climate_temp_new;                       // vector of new temperatures
  // 
  // // new ancestry data for predictions
  // int ancestry_N_new;                                           // number of new ancestry values
  // vector[ancestry_N_new] ancestry_aa_new;                       // vector of proportion aa ancestry 
  // matrix[ancestry_N_new, 2] ancestry_X_new;                     // matrix of new of continuous covariates
  // int<lower=0, upper=1> ancestry_Virus_new[ancestry_N_new];     // new categorical variable
  // 
  // // survey data
  // int surveys_N;                                                // number of survey sites
  // matrix[surveys_N, 2] surveys_X;                               // survey covariates
  // int<lower=0, upper=1> surveys_Virus[surveys_N];               // virus
  // vector[surveys_N] surveys_aa;                                 // ancestry values
  // vector[surveys_N] surveys_temp;                               // temperature values
  // 
  // // map data
  // int map_N;                                                    // number of pixels
  // matrix[map_N, 2] map_X;                                       // survey covariates
  // int<lower=0, upper=1> map_Virus[map_N];                       // virus
  // vector[map_N] map_aa;                                         // ancestry values
  // vector[map_N] map_temp;                                       // temperature values

}

parameters {
  
  // omega (probability of biting a human given a bite)
  real<lower=0> omega_ancestry_constant;                        // parameter c, lower limit
  real<lower=0, upper=1> omega_ancestry_d;                      // parameter d, upper limit
  real<lower=0> omega_ancestry_e;                               // parameter e, dose responding to halfway between c and d
  real<lower=0> omega_ancestry_sigma;                           // noise

  // alpha (biting rate)
  real<lower=0, upper=0.01> alpha_climate_constant;             // parameter c
  real<lower=0, upper=24> alpha_climate_Tmin;                   // parameter Tmin
  real<lower=24, upper=45> alpha_climate_Tmax;                  // parameter Tmax
  real<lower=0> alpha_climate_sigma;                            // noise

  // b (prob mosquito infectiousness)
  real<lower=0, upper=0.01> b_climate_constant;                 // parameter c
  real<lower=0, upper=24> b_climate_Tmin;                       // parameter Tmin
  real<lower=24, upper=45> b_climate_Tmax;                      // parameter Tmax
  real<lower=0> b_climate_sigma;                                // noise

  // // pMI denv (prob mosquito infection)
  // real<lower=0, upper=0.01> pMI_climate_constant;   // parameter c
  // real<lower=0, upper=24> pMI_climate_Tmin;         // parameter Tmin
  // real<lower=24, upper=45> pMI_climate_Tmax;        // parameter Tmax
  // real<lower=0> pMI_climate_sigma;                  // noise
  
  // pMI climate (prob mosquito infection)
  real<lower=0, upper=10> pMI_climate_rmax;                     // parameter rmax
  real<lower=10, upper=40> pMI_climate_Topt;                    // parameter Topt
  real<lower=0, upper=10> pMI_climate_a;                        // parameter a
  real<lower=0> pMI_climate_sigma;                              // noise

  // pMI ancestry (prob mosquito infection)
  real pMI_ancestry_b0;                                         // intercept
  vector[2] pMI_ancestry_beta;                                  // coefficients for continuous covariates
  real pMI_ancestry_gamma;                                      // coefficient for categorical covariate

  // EIR (extrinsic incubation rate)
  real<lower=0, upper=0.01> EIR_climate_constant;               // parameter c
  real<lower=0, upper=24> EIR_climate_Tmin;                     // parameter Tmin
  real<lower=24, upper=100> EIR_climate_Tmax;                   // parameter Tmax
  real<lower=0> EIR_climate_sigma;                              // noise

  // lifespan (1/mosquito mortality rate)
  real<lower=-0.5, upper=0.5> lf_climate_constant;              // parameter c
  real<lower=0, upper=24> lf_climate_Tmin;                      // parameter Tmin
  real<lower=24, upper=45> lf_climate_Tmax;                     // parameter Tmax
  real<lower=0> lf_climate_sigma;                               // noise
}

model {                                                         // Fit models to observed data
  // Briere equation: c*x*(x-Tmin)*sqrt(Tmax-x)
  // Quadratic equation: c*(x-Tmax)*(x-Tmin)
  // Guassian equation: rmax*exp(-0.5*(abs(x-Topt)/a)^2)
  
  // omega (probability of biting a human given a bite)
  omega_ancestry_constant ~ uniform(0,1);                       // prior for c
  omega_ancestry_d ~ uniform(0,10);                             // prior for d
  omega_ancestry_e ~ uniform(0,10);                             // prior for e
  omega_ancestry_sigma ~ uniform(0,100);                        // prior for sigma

  for(h in 1:omega_ancestry_N){
    real omega_ancestry_mu = omega_ancestry_constant + ((omega_ancestry_d - omega_ancestry_constant)/(1 + (omega_ancestry_e / omega_ancestry_aa[h])));
    omega_ancestry[h] ~ normal(omega_ancestry_mu, omega_ancestry_sigma);
  }
  
  // alpha (biting rate)
  alpha_climate_constant ~ normal(2.02E-04,0.01);               // prior for c
  alpha_climate_Tmin ~ normal(13.35,20);                        // prior for Tmin
  alpha_climate_Tmax ~ normal(40.08,20);                        // prior for Tmax
  alpha_climate_sigma ~ uniform(0,100);                         // prior for sigma

  for(i in 1:alpha_climate_N){                            
    real alpha_climate_mu = alpha_climate_constant * alpha_climate_temp[i] * (alpha_climate_temp[i] - alpha_climate_Tmin) * sqrt(alpha_climate_Tmax - alpha_climate_temp[i]);
    alpha_climate[i] ~ normal(alpha_climate_mu, alpha_climate_sigma);
  }
  if (alpha_climate_Tmin > alpha_climate_Tmax) {
    target += positive_infinity();
  }

  // b (prob mosquito infectiousness)
  b_climate_constant ~ normal(8.49E-04,0.01);                   // prior for c
  b_climate_Tmin ~ normal(17.05,20);                            // prior for Tmin
  b_climate_Tmax ~ normal(35.83,20);                            // prior for Tmax
  b_climate_sigma ~ uniform(0,100);                             // prior for sigma

  for(j in 1:b_climate_N){                          
    real b_climate_mu = b_climate_constant * b_climate_temp[j] * (b_climate_temp[j] - b_climate_Tmin) * sqrt(b_climate_Tmax - b_climate_temp[j]);
    b_climate[j] ~ normal(b_climate_mu, b_climate_sigma);
  }
  if (b_climate_Tmin > b_climate_Tmax) {
    target += positive_infinity();
  }

  // // pMI denv (prob mosquito infection)
  // pMI_climate_constant ~ normal(4.91E-04,0.01);     // prior for c
  // pMI_climate_Tmin ~ normal(12.22,20);              // prior for Tmin
  // pMI_climate_Tmax ~ normal(37.46,20);              // prior for Tmax
  // pMI_climate_sigma ~ uniform(0,100);               // prior for sigma
  // 
  // for(k in 1:pMI_climate_N){                        // Briere model
  //   real pMI_climate_mu = pMI_climate_constant * pMI_climate_temp[k] * (pMI_climate_temp[k] - pMI_climate_Tmin) * sqrt(pMI_climate_Tmax - pMI_climate_temp[k]);
  //   pMI_climate[k] ~ normal(pMI_climate_mu, pMI_climate_sigma);
  // }
  // if (pMI_climate_Tmin > pMI_climate_Tmax) {
  //   target += positive_infinity();
  // }
  
  // pMI climate (prob mosquito infection)
  pMI_climate_rmax ~ normal(0.24, 0.03);                        // prior for rmax
  pMI_climate_Topt ~ normal(30.08, 0.38);                       // prior for Topt
  pMI_climate_a ~ normal(3.60, 0.41);                           // prior for a
  pMI_climate_sigma ~ uniform(0,100);                           // prior for sigma

  for(k in 1:pMI_climate_N){                        
    real pMI_climate_mu = pMI_climate_rmax * exp(-0.5 * (fabs(pMI_climate_temp[k] - pMI_climate_Topt)/pMI_climate_a)^2);
    pMI_climate[k] ~ normal(pMI_climate_mu, pMI_climate_sigma);
  }

  // pMI ancestry (prob mosquito infection)
  pMI_ancestry_b0 ~ normal(-13,10);                             // prior for intercept
  pMI_ancestry_beta ~ normal(0, 5);                             // priors for continuous variable coefficients
  pMI_ancestry_gamma ~ normal(0, 5);                            // prior for categorical variable coefficients
  
  for (k3 in 1:pMI_ancestry_N) {
    pMI_ancestry_Infected[k3] ~ binomial(pMI_ancestry_Trials[k3], inv_logit(pMI_ancestry_b0 + pMI_ancestry_X[k3] * pMI_ancestry_beta + pMI_ancestry_gamma * pMI_ancestry_Virus[k3]));
  }

  // EIR (extrinsic incubation rate)
  EIR_climate_constant ~ normal(6.65E-05,0.01);                // prior for c
  EIR_climate_Tmin ~ normal(10.68,20);                         // prior for Tmin
  EIR_climate_Tmax ~ normal(45.90,20);                         // prior for Tmax
  EIR_climate_sigma ~ uniform(0,100);                          // prior for sigma

  for(l in 1:EIR_climate_N){                        
    real EIR_climate_mu = EIR_climate_constant * EIR_climate_temp[l] * (EIR_climate_temp[l] - EIR_climate_Tmin) * sqrt(EIR_climate_Tmax - EIR_climate_temp[l]);
    EIR_climate[l] ~ normal(EIR_climate_mu, EIR_climate_sigma);
  }
  if (EIR_climate_Tmin > EIR_climate_Tmax) {
    target += positive_infinity();
  }

  // lifespan (1/mosquito mortality rate)
  lf_climate_constant ~ normal(-1.48E-01,0.1);                 // prior for c
  lf_climate_Tmin ~ normal(9.16,20);                           // prior for Tmin
  lf_climate_Tmax ~ normal(37.73,20);                          // prior for Tmax
  lf_climate_sigma ~ uniform(0,100);                           // prior for sigma

  for(m in 1:lf_climate_N){                        
    real lf_climate_mu = lf_climate_constant * (lf_climate_temp[m] - lf_climate_Tmax) * (lf_climate_temp[m] - lf_climate_Tmin);
    lf_climate[m] ~ normal(lf_climate_mu, lf_climate_sigma);
  }
  if (lf_climate_Tmin > lf_climate_Tmax) {
    target += positive_infinity();
  }

}

generated quantities {
  
  // parameters for posterior predictive checks
  real omega_ancestry_ppc[omega_ancestry_N];                        // omega (prob biting human)
  real alpha_climate_ppc[alpha_climate_N];                          // alpha (biting rate)        
  real b_climate_ppc[b_climate_N];                                  // b (prob mosquito infectiousness)
  real pMI_climate_ppc[pMI_climate_N];                              // pMI (prob mosquito infection)
  real EIR_climate_ppc[EIR_climate_N];                              // EIR (extrinsic incubation rate)
  real lf_climate_ppc[lf_climate_N];                                // lifespan (1/mosquito mortality rate)
  real<lower=0> pMI_ancestry_Infected_ppc[pMI_ancestry_N];          //
  real<lower=0, upper=1> pMI_ancestry_propInf_ppc[pMI_ancestry_N];  //

  // new parameters for predictions
  vector[N_new] NmNh_new;                                           // ratio of mosquitoes to humans
  vector[N_new] delta_new;                                          // intrinsic incubation period
  vector[N_new] mu_h_new;                                           // human mortality rate
  vector[N_new] gamma_new;                                          // Human infectivity period
  vector[N_new] alpha_new;                                          // alpha (biting rate)
  vector[N_new] b_new;                                              // b (prob mosquito infectiousness)
  vector[N_new] EIR_new;                                            // EIR (extrinsic incubation rate)
  vector[N_new] lf_new;                                             // lifespan (1/mosquito mortality rate)
  vector[N_new] alpha_climate_new;                                  // alpha (biting rate)
  vector[N_new] b_climate_new;                                      // b (prob mosquito infectiousness)
  vector[N_new] pMI_climate_new;                                    // pMI (prob mosquito infection)
  vector[N_new] EIR_climate_new;                                    // EIR (extrinsic incubation rate)
  vector[N_new] lf_climate_new;                                     // lifespan (1/mosquito mortality rate)
  vector[N_new] pMI_ancestry_new;                                   // pMI (prob mosquito infection)
  vector[N_new] omega_new;                                          // omega (prob biting human)

  // R0 models
  vector[N_new] R0_climate_new;
  vector[N_new] R0_ancestry_new;
  vector[N_new] R0_full_new;
  
  // // parameters for models to predict across temperature 
  // vector[climate_N_new] NmNh_climate;                               // ratio of mosquitoes to humans
  // vector[climate_N_new] delta_climate;                              // intrinsic incubation period
  // vector[climate_N_new] mu_h_climate;                               // human mortality rate
  // vector[climate_N_new] gamma_climate;                              // Human infectivity period 
  // vector[climate_N_new] alpha_climate_new;                          // alpha (biting rate)
  // vector[climate_N_new] b_climate_new;                              // b (prob mosquito infectiousness)
  // vector[climate_N_new] pMI_climate_new;                            // pMI (prob mosquito infection)
  // vector[climate_N_new] EIR_climate_new;                            // EIR (extrinsic incubation rate)
  // vector[climate_N_new] lf_climate_new;                             // lifespan (1/mosquito mortality rate)
  // 
  // // parameters for models to predict across ancestry values 
  // vector[ancestry_N_new] alpha_ancestry;                            // biting rate 
  // vector[ancestry_N_new] b_ancestry;                                // probability of mosquito infectiousness 
  // vector[ancestry_N_new] lf_ancestry;                               // lifepsan
  // vector[ancestry_N_new] EIR_ancestry;                              // extrinsic incubation rate
  // vector[ancestry_N_new] pMI_ancestry_constant;                     // probability of mosquito infection 
  // vector[ancestry_N_new] NmNh_ancestry;                             // ratio of mosquitoes to humans
  // vector[ancestry_N_new] delta_ancestry;                            // intrinsic incubation period
  // vector[ancestry_N_new] mu_h_ancestry;                             // human mortality rate
  // vector[ancestry_N_new] gamma_ancestry;                            // Human infectivity period 
  // vector[ancestry_N_new] omega_ancestry_new;                        // omega (prob biting human)
  // vector[ancestry_N_new] pMI_ancestry_new;                          // pMI (prob biting human)
  // 
  // // parameters for models to predict across surveys 
  // vector[surveys_N] alpha_climate_surveys;                          // alpha (biting rate)
  // vector[surveys_N] b_climate_surveys;                              // b (prob mosquito infectiousness)
  // vector[surveys_N] pMI_climate_surveys;                            // pMI (prob mosquito infection)
  // vector[surveys_N] EIR_climate_surveys;                            // EIR (extrinsic incubation rate)
  // vector[surveys_N] lf_climate_surveys;                             // lifespan (1/mosquito mortality rate)
  // vector[surveys_N] alpha_surveys;                                  // alpha (biting rate)
  // vector[surveys_N] b_surveys;                                      // b (prob mosquito infectiousness)
  // vector[surveys_N] pMI_surveys;                                    // pMI (prob mosquito infection)
  // vector[surveys_N] EIR_surveys;                                    // EIR (extrinsic incubation rate)
  // vector[surveys_N] lf_surveys;                                     // lifespan (1/mosquito mortality rate)
  // vector[surveys_N] NmNh_surveys;                                   // ratio of mosquitoes to humans
  // vector[surveys_N] delta_surveys;                                  // intrinsic incubation period
  // vector[surveys_N] mu_h_surveys;                                   // human mortality rate
  // vector[surveys_N] gamma_surveys;                                  // Human infectivity period
  // vector[surveys_N] omega_surveys;
  // vector[surveys_N] pMI_ancestry_surveys;
  // 
  // // parameters for models to predict across map locations 
  // vector[map_N] alpha_climate_map;                          // alpha (biting rate)
  // vector[map_N] b_climate_map;                              // b (prob mosquito infectiousness)
  // vector[map_N] pMI_ancestry_map;                           // pMI (prob mosquito infection)
  // vector[map_N] EIR_climate_map;                            // EIR (extrinsic incubation rate)
  // vector[map_N] lf_climate_map;                             // lifespan (1/mosquito mortality rate)
  // vector[map_N] NmNh_map;                                   // ratio of mosquitoes to humans
  // vector[map_N] delta_map;                                  // intrinsic incubation period
  // vector[map_N] mu_h_map;                                   // human mortality rate
  // vector[map_N] gamma_map;                                  // Human infectivity period
  // vector[map_N] omega_map;
  // 
  // // R0 models to predict across gradients of temperature or ancestry
  // vector[climate_N_new] R0_climate;                                 // R0 (climate)
  // vector[ancestry_N_new] R0_ancestry_omega;                         // R0 (ancestry, omega only)
  // vector[ancestry_N_new] R0_ancestry_pMI;                           // R0 (ancestry, pMI only)
  // vector[ancestry_N_new] R0_ancestry;                               // R0 (ancestry)
  // 
  // // R0 models to predict at survey locations
  // vector[surveys_N] R0_climate_surveys;                             // R0 (climate)
  // vector[surveys_N] R0_ancestry_surveys;                            // R0 (ancestry)
  // vector[surveys_N] R0_full_surveys;                                // R0 (climate + ancestry)
  // 
  // // R0 models to predict across map locations
  // vector[map_N] R0_full_map;
  
  // posterior predictive
  // ppc omega (prob biting human)
  for (o in 1:omega_ancestry_N){
    real omega_ancestry_mu_ppc = omega_ancestry_constant + ((omega_ancestry_d - omega_ancestry_constant)/(1 + (omega_ancestry_e / omega_ancestry_aa[o])));
    omega_ancestry_ppc[o] = normal_rng(omega_ancestry_mu_ppc, omega_ancestry_sigma);
  }

  // ppc alpha (biting rate)
  for (p in 1:alpha_climate_N){                         
    real alpha_climate_mu_ppc = alpha_climate_constant * alpha_climate_temp[p] * (alpha_climate_temp[p] - alpha_climate_Tmin) * sqrt(alpha_climate_Tmax - alpha_climate_temp[p]);
    alpha_climate_ppc[p] = normal_rng(alpha_climate_mu_ppc, alpha_climate_sigma);
  }

  // ppc b (prob mosquito infectiousness)
  for (q in 1:b_climate_N){
    real b_climate_mu_ppc = b_climate_constant * b_climate_temp[q] * (b_climate_temp[q] - b_climate_Tmin) * sqrt(b_climate_Tmax - b_climate_temp[q]);
    b_climate_ppc[q] = normal_rng(b_climate_mu_ppc, b_climate_sigma);
  }

  // // ppc pMI denv (prob mosquito infection)
  // for (r in 1:pMI_climate_N){
  //   real pMI_climate_mu_ppc = pMI_climate_constant * pMI_climate_temp[r] * (pMI_climate_temp[r] - pMI_climate_Tmin) * sqrt(pMI_climate_Tmax - pMI_climate_temp[r]);
  //   pMI_climate_ppc[r] = normal_rng(pMI_climate_mu_ppc, pMI_climate_sigma);
  // }
  
 // ppc pMI climate (prob mosquito infection)
  for (r in 1:pMI_climate_N){
    real pMI_climate_mu_ppc = pMI_climate_rmax * exp(-0.5 * (fabs(pMI_climate_temp[r] - pMI_climate_Topt)/pMI_climate_a)^2);;
    pMI_climate_ppc[r] = normal_rng(pMI_climate_mu_ppc, pMI_climate_sigma);
  }

 // ppc pMI ancestry (prob mosquito infection)
  for (rr in 1:pMI_ancestry_N){
    pMI_ancestry_propInf_ppc[rr] = inv_logit(pMI_ancestry_b0 + pMI_ancestry_X[rr] * pMI_ancestry_beta + pMI_ancestry_gamma * pMI_ancestry_Virus[rr]);
    pMI_ancestry_Infected_ppc[rr] = binomial_rng(pMI_ancestry_Trials[rr], pMI_ancestry_propInf_ppc[rr]);
  }

  // ppc EIR (extrinsic incubation rate)
  for (s in 1:EIR_climate_N){                       
    real EIR_climate_mu_ppc = EIR_climate_constant * EIR_climate_temp[s] * (EIR_climate_temp[s] - EIR_climate_Tmin) * sqrt(EIR_climate_Tmax - EIR_climate_temp[s]);
    EIR_climate_ppc[s] = normal_rng(EIR_climate_mu_ppc, EIR_climate_sigma);
  }

  // ppc lifespan (1/mosquito mortality rate)
  for (u in 1:lf_climate_N){                       
    real lf_climate_mu_ppc = lf_climate_constant * (lf_climate_temp[u] - lf_climate_Tmax) * (lf_climate_temp[u] - lf_climate_Tmin);
    lf_climate_ppc[u] = normal_rng(lf_climate_mu_ppc, lf_climate_sigma);
  }

  // predict on new data
  for(zz in 1:N_new){

    // omega (prob biting | ancestry)
    omega_new[zz] = normal_rng(omega_ancestry_constant + ((omega_ancestry_d - omega_ancestry_constant)/(1 + (omega_ancestry_e / aa_new[zz]))), omega_ancestry_sigma);

    // alpha (biting rate)
    if(alpha_climate_Tmin < temp_new[zz] && alpha_climate_Tmax > temp_new[zz]){
      alpha_climate_new[zz] = normal_rng((alpha_climate_constant * temp_new[zz] * (temp_new[zz] - alpha_climate_Tmin) * sqrt(alpha_climate_Tmax - temp_new[zz])), alpha_climate_sigma);
    }
    else {
      alpha_climate_new[zz] = 0;
    }

    // b (prob mosquito infectiousness)
    if(b_climate_Tmin < temp_new[zz] && b_climate_Tmax > temp_new[zz]){
      b_climate_new[zz] = normal_rng((b_climate_constant * temp_new[zz] * (temp_new[zz] - b_climate_Tmin) * sqrt(b_climate_Tmax - temp_new[zz])), b_climate_sigma);
    }
    else {
      b_climate_new[zz] = 0;
    }

    // // pMI denv (prob mosquito infection)
    // if(pMI_climate_Tmin < temp_new[zz] && pMI_climate_Tmax > temp_new[zz]){
    //   pMI_climate_new[zz] = normal_rng((pMI_climate_constant * temp_new[zz] * (temp_new[zz] - pMI_climate_Tmin) * sqrt(pMI_climate_Tmax - temp_new[zz])), pMI_climate_sigma);
    // }
    // else {
    //   pMI_climate_new[zz] = 0;
    // }

    // pMI (prob mosquito infection)
    pMI_climate_new[zz] = normal_rng(pMI_climate_rmax * exp(-0.5 * (fabs(temp_new[zz] - pMI_climate_Topt)/pMI_climate_a)^2), pMI_climate_sigma);

    // pMI ancestry (prob mosquito infection)
    pMI_ancestry_new[zz] = inv_logit(pMI_ancestry_b0 + X_new[zz] * pMI_ancestry_beta + pMI_ancestry_gamma * Virus_new[zz]);

    // EIR (extrinsic incubation rate)
    if(EIR_climate_Tmin < temp_new[zz] && EIR_climate_Tmax > temp_new[zz]){
      EIR_climate_new[zz] = normal_rng((EIR_climate_constant * temp_new[zz] * (temp_new[zz] - EIR_climate_Tmin) * sqrt(EIR_climate_Tmax - temp_new[zz])), EIR_climate_sigma);
    }
    else {
      EIR_climate_new[zz] = 0;
    }

    // lifespan (1/mosquito mortality rate)
    if(lf_climate_Tmin < temp_new[zz] && lf_climate_Tmax > temp_new[zz]){
      lf_climate_new[zz] = normal_rng((lf_climate_constant * (temp_new[zz] - lf_climate_Tmax) * (temp_new[zz] - lf_climate_Tmin)), lf_climate_sigma);
    }
    else {
      lf_climate_new[zz] = 0;
    }

    NmNh_new[zz] = normal_rng(2,0.5);
    mu_h_new[zz] = normal_rng(4.25E-05, 0.00005);
    // delta[zz] = 1/weibull_rng(2.69, 6.70);           //zikv from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5403043/#SD1
    delta_new[zz] = 1/gamma_rng(5.9, 0.5);                  //zikv from https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0004726
    gamma_new[zz] = 1/gamma_rng(5.0, 0.5);

    alpha_new[zz] = normal_rng(0.31, 0.05);
    b_new[zz] = normal_rng(0.77, 0.05);
    EIR_new[zz] = normal_rng(0.15, 0.05);
    lf_new[zz] = normal_rng(24.36, 2);
    // pMI_new_constant[zz] = normal_rng(0.69, 0.05);


    // R0 climate
    R0_climate_new[zz] = sqrt((alpha_climate_new[zz] * b_climate_new[zz] *
    (EIR_climate_new[zz] / ((1/lf_climate_new[zz]) * ((1/lf_climate_new[zz]) + EIR_climate_new[zz])))) *
    (alpha_climate_new[zz] * pMI_climate_new[zz] * NmNh_new[zz] * (delta_new[zz] / ((delta_new[zz] + mu_h_new[zz]) *
    (gamma_new[zz] + mu_h_new[zz])))));

    // R0 ancestry (omega and pMI)
    R0_ancestry_new[zz] = sqrt(
      (omega_new[zz] * alpha_new[zz] * b_new[zz] * (EIR_new[zz] / ((1/lf_new[zz]) * ((1/lf_new[zz]) + EIR_new[zz])))) *
      (alpha_new[zz] * pMI_ancestry_new[zz] * NmNh_new[zz] * (delta_new[zz] / ((delta_new[zz] + mu_h_new[zz]) * (gamma_new[zz] + mu_h_new[zz]))))
      );

    // R0 climate + ancestry (pMI ancestry)
    R0_full_new[zz] = sqrt((omega_new[zz] * alpha_climate_new[zz] * b_climate_new[zz] *
    (EIR_climate_new[zz] / ((1/lf_climate_new[zz]) * ((1/lf_climate_new[zz]) + EIR_climate_new[zz])))) *
    (alpha_climate_new[zz] * pMI_ancestry_new[zz] * NmNh_new[zz] * (delta_new[zz] / ((delta_new[zz] + mu_h_new[zz]) *
    (gamma_new[zz] + mu_h_new[zz])))));

  }
  // 
  // // predict across new temperatures
  // for (xx in 1:climate_N_new){                      
  //   
  //   // parameters for R0 model
  //   NmNh_climate[xx] = normal_rng(2,0.5);
  //   mu_h_climate[xx] = normal_rng(4.25E-05, 0.00005);
  //   // delta[xx] = 1/weibull_rng(2.69, 6.70);                         //zikv from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5403043/#SD1
  //   delta_climate[xx] = 1/gamma_rng(5.9, 0.5);                        //zikv from https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0004726              
  //   gamma_climate[xx] = 1/gamma_rng(5.0, 0.5);
  //   
  //   // alpha (biting rate)
  //   if(alpha_climate_Tmin < climate_temp_new[xx] && alpha_climate_Tmax > climate_temp_new[xx]){
  //     alpha_climate_new[xx] = normal_rng((alpha_climate_constant * climate_temp_new[xx] * (climate_temp_new[xx] - alpha_climate_Tmin) * sqrt(alpha_climate_Tmax - climate_temp_new[xx])), alpha_climate_sigma);
  //   } 
  //   else {
  //     alpha_climate_new[xx] = 0;
  //   }
  // 
  //   // b (prob mosquito infectiousness)
  //   if(b_climate_Tmin < climate_temp_new[xx] && b_climate_Tmax > climate_temp_new[xx]){
  //     b_climate_new[xx] = normal_rng((b_climate_constant * climate_temp_new[xx] * (climate_temp_new[xx] - b_climate_Tmin) * sqrt(b_climate_Tmax - climate_temp_new[xx])), b_climate_sigma);
  //   }
  //   else {
  //     b_climate_new[xx] = 0;
  //   }
  // 
  //   // // pMI denv (prob mosquito infection)
  //   // if(pMI_climate_Tmin < climate_temp_new[xx] && pMI_climate_Tmax > climate_temp_new[xx]){
  //   //   pMI_climate_new[xx] = normal_rng((pMI_climate_constant * climate_temp_new[xx] * (climate_temp_new[xx] - pMI_climate_Tmin) * sqrt(pMI_climate_Tmax - climate_temp_new[xx])), pMI_climate_sigma);
  //   // }
  //   // else {
  //   //   pMI_climate_new[xx] = 0;
  //   // }
  //   
  //   // pMI (prob mosquito infection)
  //   pMI_climate_new[xx] = normal_rng(pMI_climate_rmax * exp(-0.5 * (fabs(climate_temp_new[xx] - pMI_climate_Topt)/pMI_climate_a)^2), pMI_climate_sigma);
  // 
  //   
  //   // EIR (extrinsic incubation rate)
  //   if(EIR_climate_Tmin < climate_temp_new[xx] && EIR_climate_Tmax > climate_temp_new[xx]){
  //     EIR_climate_new[xx] = normal_rng((EIR_climate_constant * climate_temp_new[xx] * (climate_temp_new[xx] - EIR_climate_Tmin) * sqrt(EIR_climate_Tmax - climate_temp_new[xx])), EIR_climate_sigma);
  //   } 
  //   else {
  //     EIR_climate_new[xx] = 0;
  //   }
  //   
  //   // lifespan (1/mosquito mortality rate)
  //   if(lf_climate_Tmin < climate_temp_new[xx] && lf_climate_Tmax > climate_temp_new[xx]){
  //     lf_climate_new[xx] = normal_rng((lf_climate_constant * (climate_temp_new[xx] - lf_climate_Tmax) * (climate_temp_new[xx] - lf_climate_Tmin)), lf_climate_sigma);
  //   } 
  //   else {
  //     lf_climate_new[xx] = 0;
  //   }
  //   
  //   R0_climate[xx] = sqrt((alpha_climate_new[xx] * b_climate_new[xx] *
  //   (EIR_climate_new[xx] / ((1/lf_climate_new[xx]) * ((1/lf_climate_new[xx]) + EIR_climate_new[xx])))) *
  //   (alpha_climate_new[xx] * pMI_climate_new[xx] * NmNh_climate[xx] * (delta_climate[xx] / ((delta_climate[xx] + mu_h_climate[xx]) *
  //   (gamma_climate[xx] + mu_h_climate[xx])))));
  // }
  // 
  // // predict across ancestry values
  // for (yy in 1:ancestry_N_new){
  //   omega_ancestry_new[yy] = normal_rng(omega_ancestry_constant + ((omega_ancestry_d - omega_ancestry_constant)/(1 + (omega_ancestry_e / ancestry_aa_new[yy]))), omega_ancestry_sigma);
  //   pMI_ancestry_new[yy] = inv_logit(pMI_ancestry_b0 + ancestry_X_new[yy] * pMI_ancestry_beta + pMI_ancestry_gamma * ancestry_Virus_new[yy]);
  // 
  //   alpha_ancestry[yy] = normal_rng(0.31, 0.05);
  //   b_ancestry[yy] = normal_rng(0.77, 0.05);
  //   EIR_ancestry[yy] = normal_rng(0.15, 0.05);
  //   lf_ancestry[yy] = normal_rng(24.36, 2);
  //   pMI_ancestry_constant[yy] = normal_rng(0.69, 0.05);
  //   NmNh_ancestry[yy] = normal_rng(2,0.5);
  //   mu_h_ancestry[yy] = normal_rng(4.25E-05, 0.00005);
  //   delta_ancestry[yy] = 1/gamma_rng(5.9, 0.5);                  //zikv from https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0004726              
  //   gamma_ancestry[yy] = 1/gamma_rng(5.0, 0.5);
  // 
  //  // R0 omega only
  //   R0_ancestry_omega[yy] = sqrt(
  //     (omega_ancestry_new[yy] * alpha_ancestry[yy] * b_ancestry[yy] * (EIR_ancestry[yy] / ((1/lf_ancestry[yy]) * ((1/lf_ancestry[yy]) + EIR_ancestry[yy])))) *
  //     (alpha_ancestry[yy] * pMI_ancestry_constant[yy] * NmNh_ancestry[yy] * (delta_ancestry[yy] / ((delta_ancestry[yy] + mu_h_ancestry[yy]) * (gamma_ancestry[yy] + mu_h_ancestry[yy]))))
  //     );
  //  
  //  // R0 pMI only
  //   R0_ancestry_pMI[yy] = sqrt(
  //     (alpha_ancestry[yy] * b_ancestry[yy] * (EIR_ancestry[yy] / ((1/lf_ancestry[yy]) * ((1/lf_ancestry[yy]) + EIR_ancestry[yy])))) *
  //     (alpha_ancestry[yy] * pMI_ancestry_new[yy] * NmNh_ancestry[yy] * (delta_ancestry[yy] / ((delta_ancestry[yy] + mu_h_ancestry[yy]) * (gamma_ancestry[yy] + mu_h_ancestry[yy]))))
  //     );
  //     
  //  // R0 omega and pMI
  //   R0_ancestry[yy] = sqrt(
  //     (omega_ancestry_new[yy] * alpha_ancestry[yy] * b_ancestry[yy] * (EIR_ancestry[yy] / ((1/lf_ancestry[yy]) * ((1/lf_ancestry[yy]) + EIR_ancestry[yy])))) *
  //     (alpha_ancestry[yy] * pMI_ancestry_new[yy] * NmNh_ancestry[yy] * (delta_ancestry[yy] / ((delta_ancestry[yy] + mu_h_ancestry[yy]) * (gamma_ancestry[yy] + mu_h_ancestry[yy]))))
  //     );
  // }
  // 
  // for(zz in 1:surveys_N){
  // 
  //   // omega (prob biting | ancestry)
  //   omega_surveys[zz] = normal_rng(omega_ancestry_constant + ((omega_ancestry_d - omega_ancestry_constant)/(1 + (omega_ancestry_e / surveys_aa[zz]))), omega_ancestry_sigma);
  // 
  //   // alpha (biting rate)
  //   if(alpha_climate_Tmin < surveys_temp[zz] && alpha_climate_Tmax > surveys_temp[zz]){
  //     alpha_climate_surveys[zz] = normal_rng((alpha_climate_constant * surveys_temp[zz] * (surveys_temp[zz] - alpha_climate_Tmin) * sqrt(alpha_climate_Tmax - surveys_temp[zz])), alpha_climate_sigma);
  //   }
  //   else {
  //     alpha_climate_surveys[zz] = 0;
  //   }
  // 
  //   // b (prob mosquito infectiousness)
  //   if(b_climate_Tmin < surveys_temp[zz] && b_climate_Tmax > surveys_temp[zz]){
  //     b_climate_surveys[zz] = normal_rng((b_climate_constant * surveys_temp[zz] * (surveys_temp[zz] - b_climate_Tmin) * sqrt(b_climate_Tmax - surveys_temp[zz])), b_climate_sigma);
  //   }
  //   else {
  //     b_climate_surveys[zz] = 0;
  //   }
  // 
  //   // // pMI denv (prob mosquito infection)
  //   // if(pMI_climate_Tmin < surveys_temp[zz] && pMI_climate_Tmax > surveys_temp[zz]){
  //   //   pMI_climate_new[zz] = normal_rng((pMI_climate_constant * surveys_temp[zz] * (surveys_temp[zz] - pMI_climate_Tmin) * sqrt(pMI_climate_Tmax - surveys_temp[zz])), pMI_climate_sigma);
  //   // }
  //   // else {
  //   //   pMI_climate_new[zz] = 0;
  //   // }
  // 
  //   // pMI (prob mosquito infection)
  //   pMI_climate_surveys[zz] = normal_rng(pMI_climate_rmax * exp(-0.5 * (fabs(surveys_temp[zz] - pMI_climate_Topt)/pMI_climate_a)^2), pMI_climate_sigma);
  // 
  //   // pMI ancestry (prob mosquito infection)
  //   pMI_ancestry_surveys[zz] = inv_logit(pMI_ancestry_b0 + surveys_X[zz] * pMI_ancestry_beta + pMI_ancestry_gamma * surveys_Virus[zz]);
  // 
  //   // EIR (extrinsic incubation rate)
  //   if(EIR_climate_Tmin < surveys_temp[zz] && EIR_climate_Tmax > surveys_temp[zz]){
  //     EIR_climate_surveys[zz] = normal_rng((EIR_climate_constant * surveys_temp[zz] * (surveys_temp[zz] - EIR_climate_Tmin) * sqrt(EIR_climate_Tmax - surveys_temp[zz])), EIR_climate_sigma);
  //   }
  //   else {
  //     EIR_climate_surveys[zz] = 0;
  //   }
  // 
  //   // lifespan (1/mosquito mortality rate)
  //   if(lf_climate_Tmin < surveys_temp[zz] && lf_climate_Tmax > surveys_temp[zz]){
  //     lf_climate_surveys[zz] = normal_rng((lf_climate_constant * (surveys_temp[zz] - lf_climate_Tmax) * (surveys_temp[zz] - lf_climate_Tmin)), lf_climate_sigma);
  //   }
  //   else {
  //     lf_climate_surveys[zz] = 0;
  //   }
  // 
  //   NmNh_surveys[zz] = normal_rng(2,0.5);
  //   mu_h_surveys[zz] = normal_rng(4.25E-05, 0.00005);
  //   // delta[zz] = 1/weibull_rng(2.69, 6.70);           //zikv from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5403043/#SD1
  //   delta_surveys[zz] = 1/gamma_rng(5.9, 0.5);                  //zikv from https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0004726
  //   gamma_surveys[zz] = 1/gamma_rng(5.0, 0.5);
  // 
  //   alpha_surveys[zz] = normal_rng(0.31, 0.05);
  //   b_surveys[zz] = normal_rng(0.77, 0.05);
  //   EIR_surveys[zz] = normal_rng(0.15, 0.05);
  //   lf_surveys[zz] = normal_rng(24.36, 2); 
  //   // pMI_surveys_constant[zz] = normal_rng(0.69, 0.05);
  // 
  // 
  //   // R0 climate
  //   R0_climate_surveys[zz] = sqrt((alpha_climate_surveys[zz] * b_climate_surveys[zz] *
  //   (EIR_climate_surveys[zz] / ((1/lf_climate_surveys[zz]) * ((1/lf_climate_surveys[zz]) + EIR_climate_surveys[zz])))) *
  //   (alpha_climate_surveys[zz] * pMI_climate_surveys[zz] * NmNh_surveys[zz] * (delta_surveys[zz] / ((delta_surveys[zz] + mu_h_surveys[zz]) *
  //   (gamma_surveys[zz] + mu_h_surveys[zz])))));
  // 
  //   // R0 ancestry (omega and pMI)
  //   R0_ancestry_surveys[zz] = sqrt(
  //     (omega_surveys[zz] * alpha_surveys[zz] * b_surveys[zz] * (EIR_surveys[zz] / ((1/lf_surveys[zz]) * ((1/lf_surveys[zz]) + EIR_surveys[zz])))) *
  //     (alpha_surveys[zz] * pMI_ancestry_surveys[zz] * NmNh_surveys[zz] * (delta_surveys[zz] / ((delta_surveys[zz] + mu_h_surveys[zz]) * (gamma_surveys[zz] + mu_h_surveys[zz]))))
  //     );
  // 
  //   // R0 climate + ancestry (pMI ancestry)
  //   R0_full_surveys[zz] = sqrt((omega_surveys[zz] * alpha_climate_surveys[zz] * b_climate_surveys[zz] *
  //   (EIR_climate_surveys[zz] / ((1/lf_climate_surveys[zz]) * ((1/lf_climate_surveys[zz]) + EIR_climate_surveys[zz])))) *
  //   (alpha_climate_surveys[zz] * pMI_ancestry_surveys[zz] * NmNh_surveys[zz] * (delta_surveys[zz] / ((delta_surveys[zz] + mu_h_surveys[zz]) *
  //   (gamma_surveys[zz] + mu_h_surveys[zz])))));
  // 
  // }
  // 
  // for(mm in 1:map_N){
  // 
  //   // omega (prob biting | ancestry)
  //   omega_map[mm] = normal_rng(omega_ancestry_constant + ((omega_ancestry_d - omega_ancestry_constant)/(1 + (omega_ancestry_e / map_aa[mm]))), omega_ancestry_sigma);
  // 
  //   // alpha (biting rate)
  //   if(alpha_climate_Tmin < map_temp[mm] && alpha_climate_Tmax > map_temp[mm]){
  //     alpha_climate_map[mm] = normal_rng((alpha_climate_constant * map_temp[mm] * (map_temp[mm] - alpha_climate_Tmin) * sqrt(alpha_climate_Tmax - map_temp[mm])), alpha_climate_sigma);
  //   }
  //   else {
  //     alpha_climate_map[mm] = 0;
  //   }
  // 
  //   // b (prob mosquito infectiousness)
  //   if(b_climate_Tmin < map_temp[mm] && b_climate_Tmax > map_temp[mm]){
  //     b_climate_map[mm] = normal_rng((b_climate_constant * map_temp[mm] * (map_temp[mm] - b_climate_Tmin) * sqrt(b_climate_Tmax - map_temp[mm])), b_climate_sigma);
  //   }
  //   else {
  //     b_climate_map[mm] = 0;
  //   }
  // 
  //   // pMI ancestry (prob mosquito infection)
  //   pMI_ancestry_map[mm] = inv_logit(pMI_ancestry_b0 + map_X[mm] * pMI_ancestry_beta + pMI_ancestry_gamma * map_Virus[mm]);
  // 
  //   // EIR (extrinsic incubation rate)
  //   if(EIR_climate_Tmin < map_temp[mm] && EIR_climate_Tmax > map_temp[mm]){
  //     EIR_climate_map[mm] = normal_rng((EIR_climate_constant * map_temp[mm] * (map_temp[mm] - EIR_climate_Tmin) * sqrt(EIR_climate_Tmax - map_temp[mm])), EIR_climate_sigma);
  //   }
  //   else {
  //     EIR_climate_map[mm] = 0;
  //   }
  // 
  //   // lifespan (1/mosquito mortality rate)
  //   if(lf_climate_Tmin < map_temp[mm] && lf_climate_Tmax > map_temp[mm]){
  //     lf_climate_map[mm] = normal_rng((lf_climate_constant * (map_temp[mm] - lf_climate_Tmax) * (map_temp[mm] - lf_climate_Tmin)), lf_climate_sigma);
  //   }
  //   else {
  //     lf_climate_map[mm] = 0;
  //   }
  // 
  //   NmNh_map[mm] = normal_rng(2,0.5);
  //   mu_h_map[mm] = normal_rng(4.25E-05, 0.00005);
  //   // delta[mm] = 1/weibull_rng(2.69, 6.70);           //zikv from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5403043/#SD1
  //   delta_map[mm] = 1/gamma_rng(5.9, 0.5);                  //zikv from https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0004726
  //   gamma_map[mm] = 1/gamma_rng(5.0, 0.5);
  // 
  //   // R0 climate + ancestry (pMI ancestry)
  //   R0_full_map[mm] = sqrt((omega_map[mm] * alpha_climate_map[mm] * b_climate_map[mm] *
  //   (EIR_climate_map[mm] / ((1/lf_climate_map[mm]) * ((1/lf_climate_map[mm]) + EIR_climate_map[mm])))) *
  //   (alpha_climate_map[mm] * pMI_ancestry_map[mm] * NmNh_map[mm] * (delta_map[mm] / ((delta_map[mm] + mu_h_map[mm]) *
  //   (gamma_map[mm] + mu_h_map[mm])))));
  // 
  // }

}
