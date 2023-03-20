data {
  // alpha (biting rate)
  int alpha_climate_N;                              // number of observations
  vector[alpha_climate_N] alpha_climate_temp;       // vector of temperatures
  vector[alpha_climate_N] alpha_climate;            // vector of trait b 

  // b denv (prob mosquito infectiousness)
  int b_denv_climate_N;                             // number of observations
  vector[b_denv_climate_N] b_denv_climate_temp;     // vector of temperatures
  vector[b_denv_climate_N] b_denv_climate;          // vector of trait b 
  
  // b zikv (prob mosquito infectiousness)
  int b_zikv_climate_N;                             // number of observations
  vector[b_zikv_climate_N] b_zikv_climate_temp;     // vector of temperatures
  vector[b_zikv_climate_N] b_zikv_climate;          // vector of trait b 

  // pMI denv (prob mosquito infection)
  int pMI_denv_climate_N;                           // number of observations
  vector[pMI_denv_climate_N] pMI_denv_climate_temp; // vector of temperatures
  vector[pMI_denv_climate_N] pMI_denv_climate;      // vector of trait pMI

  // pMI zikv (prob mosquito infection)
  int pMI_zikv_climate_N;                           // number of observations
  vector[pMI_zikv_climate_N] pMI_zikv_climate_temp; // vector of temperatures
  vector[pMI_zikv_climate_N] pMI_zikv_climate;      // vector of trait pMI

  // EIR (extrinsic incubation rate)
  int EIR_climate_N;                                // number of observations
  vector[EIR_climate_N] EIR_climate_temp;           // vector of temperatures
  vector[EIR_climate_N] EIR_climate;                // vector of trait EIR

  // lifespan (1/mosquito mortality rate)
  int lf_climate_N;                                 // number of observations
  vector[lf_climate_N] lf_climate_temp;             // vector of temperatures
  vector[lf_climate_N] lf_climate;                  // vector of trait lf

  // new climate data for predictions 
  int climate_N_new;                                // number of new temperatures
  vector[climate_N_new] climate_temp_new;           // vector of new temperatures
}

parameters {
  // alpha (biting rate)
  real<lower=0, upper=0.01> alpha_climate_constant; // parameter c
  real<lower=0, upper=24> alpha_climate_Tmin;       // parameter Tmin
  real<lower=24, upper=45> alpha_climate_Tmax;      // parameter Tmax
  real<lower=0> alpha_climate_sigma;                // noise

  // b denv (prob mosquito infectiousness)
  real<lower=0, upper=0.01> b_denv_climate_constant;// parameter c
  real<lower=0, upper=24> b_denv_climate_Tmin;      // parameter Tmin
  real<lower=24, upper=45> b_denv_climate_Tmax;     // parameter Tmax
  real<lower=0> b_denv_climate_sigma;               // noise

  // b zikv (prob mosquito infectiousness)
  real<lower=0, upper=0.01> b_zikv_climate_constant;// parameter c
  real<lower=0, upper=24> b_zikv_climate_Tmin;      // parameter Tmin
  real<lower=24, upper=45> b_zikv_climate_Tmax;     // parameter Tmax
  real<lower=0> b_zikv_climate_sigma;               // noise

  // pMI denv (prob mosquito infection)
  real<lower=0, upper=0.01> pMI_denv_climate_constant;// parameter c
  real<lower=0, upper=24> pMI_denv_climate_Tmin;    // parameter Tmin
  real<lower=24, upper=45> pMI_denv_climate_Tmax;   // parameter Tmax
  real<lower=0> pMI_denv_climate_sigma;             // noise

  // pMI zikv (prob mosquito infection)
  real<lower=0, upper=0.01> pMI_zikv_climate_constant;// parameter c
  real<lower=0, upper=24> pMI_zikv_climate_Tmin;    // parameter Tmin
  real<lower=24, upper=45> pMI_zikv_climate_Tmax;   // parameter Tmax
  real<lower=0> pMI_zikv_climate_sigma;             // noise

  // EIR (extrinsic incubation rate)
  real<lower=0, upper=0.01> EIR_climate_constant;   // parameter c
  real<lower=0, upper=24> EIR_climate_Tmin;         // parameter Tmin
  real<lower=24, upper=100> EIR_climate_Tmax;       // parameter Tmax
  real<lower=0> EIR_climate_sigma;                  // noise

  // lifespan (1/mosquito mortality rate)
  real<lower=-0.5, upper=0.5> lf_climate_constant;  // parameter c
  real<lower=0, upper=24> lf_climate_Tmin;          // parameter Tmin
  real<lower=24, upper=45> lf_climate_Tmax;         // parameter Tmax
  real<lower=0> lf_climate_sigma;                   // noise
}

model {                                             // Fit models to observed data
  // Briere equation: c*x*(x-Tmin)*sqrt(Tmax-x)
  // Quadratic equation: c*(x-Tmax)*(x-Tmin)

  // alpha (biting rate)
  alpha_climate_constant ~ normal(2.02E-04,0.01);   // prior for c
  alpha_climate_Tmin ~ normal(13.35,20);            // prior for Tmin
  alpha_climate_Tmax ~ normal(40.08,20);            // prior for Tmax
  alpha_climate_sigma ~ uniform(0,100);             // prior for sigma

  for(i in 1:alpha_climate_N){                      // Briere model
    real alpha_climate_mu = alpha_climate_constant * alpha_climate_temp[i] * (alpha_climate_temp[i] - alpha_climate_Tmin) * sqrt(alpha_climate_Tmax - alpha_climate_temp[i]);
    alpha_climate[i] ~ normal(alpha_climate_mu, alpha_climate_sigma);
  }
  if (alpha_climate_Tmin > alpha_climate_Tmax) {
    target += positive_infinity();
  }

  // b denv (prob mosquito infectiousness)
  b_denv_climate_constant ~ normal(8.49E-04,0.01);  // prior for c
  b_denv_climate_Tmin ~ normal(17.05,20);           // prior for Tmin
  b_denv_climate_Tmax ~ normal(35.83,20);           // prior for Tmax
  b_denv_climate_sigma ~ uniform(0,100);            // prior for sigma

  for(j in 1:b_denv_climate_N){                     // Briere model
    real b_denv_climate_mu = b_denv_climate_constant * b_denv_climate_temp[j] * (b_denv_climate_temp[j] - b_denv_climate_Tmin) * sqrt(b_denv_climate_Tmax - b_denv_climate_temp[j]);
    b_denv_climate[j] ~ normal(b_denv_climate_mu, b_denv_climate_sigma);
  }
  if (b_denv_climate_Tmin > b_denv_climate_Tmax) {
    target += positive_infinity();
  }

  // NEED TO CHANGE THESE PRIORS!
  // b zikv (prob mosquito infectiousness)
  b_zikv_climate_constant ~ normal(8.49E-04,0.01);  // prior for c
  b_zikv_climate_Tmin ~ normal(17.05,20);           // prior for Tmin
  b_zikv_climate_Tmax ~ normal(35.83,20);           // prior for Tmax
  b_zikv_climate_sigma ~ uniform(0,100);            // prior for sigma

  for(jj in 1:b_zikv_climate_N){                     // Briere model
    real b_zikv_climate_mu = b_zikv_climate_constant * b_zikv_climate_temp[jj] * (b_zikv_climate_temp[jj] - b_zikv_climate_Tmin) * sqrt(b_zikv_climate_Tmax - b_zikv_climate_temp[jj]);
    b_zikv_climate[jj] ~ normal(b_zikv_climate_mu, b_zikv_climate_sigma);
  }
  if (b_zikv_climate_Tmin > b_zikv_climate_Tmax) {
    target += positive_infinity();
  }

  
  // pMI denv (prob mosquito infection)
  pMI_denv_climate_constant ~ normal(4.91E-04,0.01);// prior for c
  pMI_denv_climate_Tmin ~ normal(12.22,20);         // prior for Tmin
  pMI_denv_climate_Tmax ~ normal(37.46,20);         // prior for Tmax
  pMI_denv_climate_sigma ~ uniform(0,100);          // prior for sigma

  for(k in 1:pMI_denv_climate_N){                   // Briere model
    real pMI_denv_climate_mu = pMI_denv_climate_constant * pMI_denv_climate_temp[k] * (pMI_denv_climate_temp[k] - pMI_denv_climate_Tmin) * sqrt(pMI_denv_climate_Tmax - pMI_denv_climate_temp[k]);
    pMI_denv_climate[k] ~ normal(pMI_denv_climate_mu, pMI_denv_climate_sigma);
  }
  if (pMI_denv_climate_Tmin > pMI_denv_climate_Tmax) {
    target += positive_infinity();
  }
  
  // NEED TO CHANGE THESE PRIORS!
  // pMI zikv (prob mosquito infection)
  pMI_zikv_climate_constant ~ normal(4.91E-04,0.01);// prior for c
  pMI_zikv_climate_Tmin ~ normal(12.22,20);         // prior for Tmin
  pMI_zikv_climate_Tmax ~ normal(37.46,20);         // prior for Tmax
  pMI_zikv_climate_sigma ~ uniform(0,100);          // prior for sigma

  for(kk in 1:pMI_zikv_climate_N){                   // Briere model
    real pMI_zikv_climate_mu = pMI_zikv_climate_constant * pMI_zikv_climate_temp[kk] * (pMI_zikv_climate_temp[kk] - pMI_zikv_climate_Tmin) * sqrt(pMI_zikv_climate_Tmax - pMI_zikv_climate_temp[kk]);
    pMI_zikv_climate[kk] ~ normal(pMI_zikv_climate_mu, pMI_zikv_climate_sigma);
  }
  if (pMI_zikv_climate_Tmin > pMI_zikv_climate_Tmax) {
    target += positive_infinity();
  }

  // EIR (extrinsic incubation rate)
  EIR_climate_constant ~ normal(6.65E-05,0.01);     // prior for c
  EIR_climate_Tmin ~ normal(10.68,20);              // prior for Tmin
  EIR_climate_Tmax ~ normal(45.90,20);              // prior for Tmax
  EIR_climate_sigma ~ uniform(0,100);               // prior for sigma

  for(l in 1:EIR_climate_N){                        // Briere model
    real EIR_climate_mu = EIR_climate_constant * EIR_climate_temp[l] * (EIR_climate_temp[l] - EIR_climate_Tmin) * sqrt(EIR_climate_Tmax - EIR_climate_temp[l]);
    EIR_climate[l] ~ normal(EIR_climate_mu, EIR_climate_sigma);
  }
  if (EIR_climate_Tmin > EIR_climate_Tmax) {
    target += positive_infinity();
  }

  // lifespan (1/mosquito mortality rate)
  lf_climate_constant ~ normal(-1.48E-01,0.1);     // prior for c
  lf_climate_Tmin ~ normal(9.16,20);               // prior for Tmin
  lf_climate_Tmax ~ normal(37.73,20);              // prior for Tmax
  lf_climate_sigma ~ uniform(0,100);               // prior for sigma

  for(m in 1:lf_climate_N){                        // Quadratic model
    real lf_climate_mu = lf_climate_constant * (lf_climate_temp[m] - lf_climate_Tmax) * (lf_climate_temp[m] - lf_climate_Tmin);
    lf_climate[m] ~ normal(lf_climate_mu, lf_climate_sigma);
  }
  if (lf_climate_Tmin > lf_climate_Tmax) {
    target += positive_infinity();
  }

}

generated quantities {
  
  // posterior predictive check values
  real alpha_climate_ppc[alpha_climate_N];          // alpha (biting rate)        
  real b_denv_climate_ppc[b_denv_climate_N];        // b (prob mosquito infectiousness)
  real b_zikv_climate_ppc[b_zikv_climate_N];        // b (prob mosquito infectiousness)
  real pMI_denv_climate_ppc[pMI_denv_climate_N];    // pMI (prob mosquito infection)
  real pMI_zikv_climate_ppc[pMI_zikv_climate_N];    // pMI (prob mosquito infection)
  real EIR_climate_ppc[EIR_climate_N];              // EIR (extrinsic incubation rate)
  real lf_climate_ppc[lf_climate_N];                // lifespan (1/mosquito mortality rate)

  // new vectors for parameters
  vector[climate_N_new] NmNh;                       // ratio of mosquitoes to humans
  vector[climate_N_new] mu_h;                       // human mortality rate
  vector[climate_N_new] gamma;                      // Human infectivity period 
  
  vector[climate_N_new] alpha_climate_new;          // alpha (biting rate)
  vector[climate_N_new] b_denv_climate_new;         // b (prob mosquito infectiousness)
  vector[climate_N_new] b_zikv_climate_new;         // b (prob mosquito infectiousness)
  vector[climate_N_new] pMI_denv_climate_new;       // pMI (prob mosquito infection)
  vector[climate_N_new] pMI_zikv_climate_new;       // pMI (prob mosquito infection)
  vector[climate_N_new] EIR_climate_new;            // EIR (extrinsic incubation rate)
  vector[climate_N_new] lf_climate_new;             // lifespand (1/mosquito mortality rate)
  vector[climate_N_new] R0_denv_climate;            // R0
  vector[climate_N_new] R0_zikv_climate;            // R0
  
  // ppc alpha (biting rate)
  for (p in 1:alpha_climate_N){                         
    real alpha_climate_mu_ppc = alpha_climate_constant * alpha_climate_temp[p] * (alpha_climate_temp[p] - alpha_climate_Tmin) * sqrt(alpha_climate_Tmax - alpha_climate_temp[p]);
    alpha_climate_ppc[p] = normal_rng(alpha_climate_mu_ppc, alpha_climate_sigma);
  }

  // ppc b denv (prob mosquito infectiousness)
  for (q in 1:b_denv_climate_N){                         
    real b_denv_climate_mu_ppc = b_denv_climate_constant * b_denv_climate_temp[q] * (b_denv_climate_temp[q] - b_denv_climate_Tmin) * sqrt(b_denv_climate_Tmax - b_denv_climate_temp[q]);
    b_denv_climate_ppc[q] = normal_rng(b_denv_climate_mu_ppc, b_denv_climate_sigma);
  }

  // ppc b zikv (prob mosqquito infectiousness)
  for (qq in 1:b_zikv_climate_N){                         
    real b_zikv_climate_mu_ppc = b_zikv_climate_constant * b_zikv_climate_temp[qq] * (b_zikv_climate_temp[qq] - b_zikv_climate_Tmin) * sqrt(b_zikv_climate_Tmax - b_zikv_climate_temp[qq]);
    b_zikv_climate_ppc[qq] = normal_rng(b_zikv_climate_mu_ppc, b_zikv_climate_sigma);
  }

  // ppc pMI denv (prob mosquito infection)
  for (r in 1:pMI_denv_climate_N){                       
    real pMI_denv_climate_mu_ppc = pMI_denv_climate_constant * pMI_denv_climate_temp[r] * (pMI_denv_climate_temp[r] - pMI_denv_climate_Tmin) * sqrt(pMI_denv_climate_Tmax - pMI_denv_climate_temp[r]);
    pMI_denv_climate_ppc[r] = normal_rng(pMI_denv_climate_mu_ppc, pMI_denv_climate_sigma);
  }

  // ppc pMI zikv (prrob mosquito infection)
  for (rr in 1:pMI_zikv_climate_N){                       
    real pMI_zikv_climate_mu_ppc = pMI_zikv_climate_constant * pMI_zikv_climate_temp[rr] * (pMI_zikv_climate_temp[rr] - pMI_zikv_climate_Tmin) * sqrt(pMI_zikv_climate_Tmax - pMI_zikv_climate_temp[rr]);
    pMI_zikv_climate_ppc[rr] = normal_rng(pMI_zikv_climate_mu_ppc, pMI_zikv_climate_sigma);
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
  
  // predict across new temperatures
  for (xx in 1:climate_N_new){                      
  
    // alpha (biting rate)
    if(alpha_climate_Tmin < climate_temp_new[xx] && alpha_climate_Tmax > climate_temp_new[xx]){
      alpha_climate_new[xx] = normal_rng((alpha_climate_constant * climate_temp_new[xx] * (climate_temp_new[xx] - alpha_climate_Tmin) * sqrt(alpha_climate_Tmax - climate_temp_new[xx])), alpha_climate_sigma);
    } 
    else {
      alpha_climate_new[xx] = 0;
    }

    // b denv (prob mosquito infectiousness)
    if(b_denv_climate_Tmin < climate_temp_new[xx] && b_denv_climate_Tmax > climate_temp_new[xx]){
      b_denv_climate_new[xx] = normal_rng((b_denv_climate_constant * climate_temp_new[xx] * (climate_temp_new[xx] - b_denv_climate_Tmin) * sqrt(b_denv_climate_Tmax - climate_temp_new[xx])), b_denv_climate_sigma);
    } 
    else {
      b_denv_climate_new[xx] = 0;
    }
    
    // b zikv (prob mosquito infectiousness)
    if(b_zikv_climate_Tmin < climate_temp_new[xx] && b_zikv_climate_Tmax > climate_temp_new[xx]){
      b_zikv_climate_new[xx] = normal_rng((b_zikv_climate_constant * climate_temp_new[xx] * (climate_temp_new[xx] - b_zikv_climate_Tmin) * sqrt(b_zikv_climate_Tmax - climate_temp_new[xx])), b_zikv_climate_sigma);
    } 
    else {
      b_zikv_climate_new[xx] = 0;
    }
    
    // pMI denv (prob mosquito infection)
    if(pMI_denv_climate_Tmin < climate_temp_new[xx] && pMI_denv_climate_Tmax > climate_temp_new[xx]){
      pMI_denv_climate_new[xx] = normal_rng((pMI_denv_climate_constant * climate_temp_new[xx] * (climate_temp_new[xx] - pMI_denv_climate_Tmin) * sqrt(pMI_denv_climate_Tmax - climate_temp_new[xx])), pMI_denv_climate_sigma);
    } 
    else {
      pMI_denv_climate_new[xx] = 0;
    }
    
    // pMI zikv (prob mosquito infection)
    if(pMI_zikv_climate_Tmin < climate_temp_new[xx] && pMI_zikv_climate_Tmax > climate_temp_new[xx]){
      pMI_zikv_climate_new[xx] = normal_rng((pMI_zikv_climate_constant * climate_temp_new[xx] * (climate_temp_new[xx] - pMI_zikv_climate_Tmin) * sqrt(pMI_zikv_climate_Tmax - climate_temp_new[xx])), pMI_zikv_climate_sigma);
    } 
    else {
      pMI_zikv_climate_new[xx] = 0;
    }
    
    // EIR (extrinsic incubation rate)
    if(EIR_climate_Tmin < climate_temp_new[xx] && EIR_climate_Tmax > climate_temp_new[xx]){
      EIR_climate_new[xx] = normal_rng((EIR_climate_constant * climate_temp_new[xx] * (climate_temp_new[xx] - EIR_climate_Tmin) * sqrt(EIR_climate_Tmax - climate_temp_new[xx])), EIR_climate_sigma);
    } 
    else {
      EIR_climate_new[xx] = 0;
    }
    
    // lifespan (1/mosquito mortality rate)
    if(lf_climate_Tmin < climate_temp_new[xx] && lf_climate_Tmax > climate_temp_new[xx]){
      lf_climate_new[xx] = normal_rng((lf_climate_constant * (climate_temp_new[xx] - lf_climate_Tmax) * (climate_temp_new[xx] - lf_climate_Tmin)), lf_climate_sigma);
    } 
    else {
      lf_climate_new[xx] = 0;
    }
    
    NmNh[xx] = normal_rng(2,0.5);
    mu_h[xx] = normal_rng(4.25E-05, 5E-06);
    gamma[xx] = normal_rng(0.2, 0.06);

    R0_denv_climate[xx] = sqrt((alpha_climate_new[xx] * b_denv_climate_new[xx] *
    (EIR_climate_new[xx] / ((1/lf_climate_new[xx]) * ((1/lf_climate_new[xx]) + EIR_climate_new[xx])))) *
    (alpha_climate_new[xx] * pMI_denv_climate_new[xx] * NmNh[xx] * (gamma[xx] / ((gamma[xx] + mu_h[xx]) *
    (gamma[xx] + mu_h[xx])))));
    
    R0_zikv_climate[xx] = sqrt((alpha_climate_new[xx] * b_zikv_climate_new[xx] *
    (EIR_climate_new[xx] / ((1/lf_climate_new[xx]) * ((1/lf_climate_new[xx]) + EIR_climate_new[xx])))) *
    (alpha_climate_new[xx] * pMI_zikv_climate_new[xx] * NmNh[xx] * (gamma[xx] / ((gamma[xx] + mu_h[xx]) *
    (gamma[xx] + mu_h[xx])))));

  }
}
