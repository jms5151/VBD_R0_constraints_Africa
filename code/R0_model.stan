data {
  // alpha (biting rate)
  int alpha_climate_N;                              // number of observations
  vector[alpha_climate_N] alpha_climate_temp;       // vector of temperatures
  vector[alpha_climate_N] alpha_climate;            // vector of trait alpha

  // b (prob mosquito infectiousness)
  int b_climate_N;                                  // number of observations
  vector[b_climate_N] b_climate_temp;               // vector of temperatures
  vector[b_climate_N] b_climate;                    // vector of trait b
  
  // pMI (prob mosquito infection)
  int pMI_climate_N;                                // number of observations
  vector[pMI_climate_N] pMI_climate_temp;           // vector of temperatures
  vector[pMI_climate_N] pMI_climate;                // vector of trait pMI

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

  // b (prob mosquito infectiousness)
  real<lower=0, upper=0.01> b_climate_constant;     // parameter c
  real<lower=0, upper=24> b_climate_Tmin;           // parameter Tmin
  real<lower=24, upper=45> b_climate_Tmax;          // parameter Tmax
  real<lower=0> b_climate_sigma;                    // noise

  // // pMI denv (prob mosquito infection)
  // real<lower=0, upper=0.01> pMI_climate_constant;   // parameter c
  // real<lower=0, upper=24> pMI_climate_Tmin;         // parameter Tmin
  // real<lower=24, upper=45> pMI_climate_Tmax;        // parameter Tmax
  // real<lower=0> pMI_climate_sigma;                  // noise
  
  // pMI (prob mosquito infection)
  real<lower=0, upper=10> pMI_climate_rmax;         // parameter rmax
  real<lower=10, upper=40> pMI_climate_Topt;        // parameter Topt
  real<lower=0, upper=10> pMI_climate_a;            // parameter a
  real<lower=0> pMI_climate_sigma;                  // noise

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
  // Guassian equation: rmax*exp(-0.5*(abs(x-Topt)/a)^2)

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

  // b (prob mosquito infectiousness)
  b_climate_constant ~ normal(8.49E-04,0.01);       // prior for c
  b_climate_Tmin ~ normal(17.05,20);                // prior for Tmin
  b_climate_Tmax ~ normal(35.83,20);                // prior for Tmax
  b_climate_sigma ~ uniform(0,100);                 // prior for sigma

  for(j in 1:b_climate_N){                          // Briere model
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
  
  // pMI (prob mosquito infection)
  pMI_climate_rmax ~ normal(0.24, 0.03);            // prior for rmax
  pMI_climate_Topt ~ normal(30.08, 0.38);           // prior for Topt
  pMI_climate_a ~ normal(3.60, 0.41);               // prior for a
  pMI_climate_sigma ~ uniform(0,100);               // prior for sigma

  for(k in 1:pMI_climate_N){                        // gaussian model
    real pMI_climate_mu = pMI_climate_rmax * exp(-0.5 * (fabs(pMI_climate_temp[k] - pMI_climate_Topt)/pMI_climate_a)^2);
    pMI_climate[k] ~ normal(pMI_climate_mu, pMI_climate_sigma);
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
  real b_climate_ppc[b_climate_N];                  // b (prob mosquito infectiousness)
  real pMI_climate_ppc[pMI_climate_N];              // pMI (prob mosquito infection)
  real EIR_climate_ppc[EIR_climate_N];              // EIR (extrinsic incubation rate)
  real lf_climate_ppc[lf_climate_N];                // lifespan (1/mosquito mortality rate)

  // new vectors for parameters
  vector[climate_N_new] NmNh;                       // ratio of mosquitoes to humans
  vector[climate_N_new] delta;                      // intrinsic incubation period
  vector[climate_N_new] mu_h;                       // human mortality rate
  vector[climate_N_new] gamma;                      // Human infectivity period 
  
  vector[climate_N_new] alpha_climate_new;          // alpha (biting rate)
  vector[climate_N_new] b_climate_new;              // b (prob mosquito infectiousness)
  vector[climate_N_new] pMI_climate_new;            // pMI (prob mosquito infection)
  vector[climate_N_new] EIR_climate_new;            // EIR (extrinsic incubation rate)
  vector[climate_N_new] lf_climate_new;             // lifespand (1/mosquito mortality rate)
  vector[climate_N_new] R0;                         // R0
  
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
  
 // ppc pMI (prob mosquito infection)
  for (r in 1:pMI_climate_N){
    real pMI_climate_mu_ppc = pMI_climate_rmax * exp(-0.5 * (fabs(pMI_climate_temp[r] - pMI_climate_Topt)/pMI_climate_a)^2);;
    pMI_climate_ppc[r] = normal_rng(pMI_climate_mu_ppc, pMI_climate_sigma);
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

    // b (prob mosquito infectiousness)
    if(b_climate_Tmin < climate_temp_new[xx] && b_climate_Tmax > climate_temp_new[xx]){
      b_climate_new[xx] = normal_rng((b_climate_constant * climate_temp_new[xx] * (climate_temp_new[xx] - b_climate_Tmin) * sqrt(b_climate_Tmax - climate_temp_new[xx])), b_climate_sigma);
    }
    else {
      b_climate_new[xx] = 0;
    }

    // // pMI denv (prob mosquito infection)
    // if(pMI_climate_Tmin < climate_temp_new[xx] && pMI_climate_Tmax > climate_temp_new[xx]){
    //   pMI_climate_new[xx] = normal_rng((pMI_climate_constant * climate_temp_new[xx] * (climate_temp_new[xx] - pMI_climate_Tmin) * sqrt(pMI_climate_Tmax - climate_temp_new[xx])), pMI_climate_sigma);
    // }
    // else {
    //   pMI_climate_new[xx] = 0;
    // }
    
    // pMI (prob mosquito infection)
    pMI_climate_new[xx] = normal_rng(pMI_climate_rmax * exp(-0.5 * (fabs(climate_temp_new[xx] - pMI_climate_Topt)/pMI_climate_a)^2), pMI_climate_sigma);

    
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
    mu_h[xx] = normal_rng(4.25E-05, 0.00005);
    // delta[xx] = 1/weibull_rng(2.69, 6.70);           //zikv from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5403043/#SD1
    delta[xx] = 1/gamma_rng(5.9, 0.5);                  //zikv from https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0004726              
    gamma[xx] = 1/gamma_rng(5.0, 0.5);

    R0[xx] = sqrt((alpha_climate_new[xx] * b_climate_new[xx] *
    (EIR_climate_new[xx] / ((1/lf_climate_new[xx]) * ((1/lf_climate_new[xx]) + EIR_climate_new[xx])))) *
    (alpha_climate_new[xx] * pMI_climate_new[xx] * NmNh[xx] * (delta[xx] / ((delta[xx] + mu_h[xx]) *
    (gamma[xx] + mu_h[xx])))));
  }
}