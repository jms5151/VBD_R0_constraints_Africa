data {

  
  // b (prob mosquito infectiousness)
  int b_climate_N;                                  // number of observations
  vector[b_climate_N] b_climate_temp;               // vector of temperatures
  vector[b_climate_N] b_climate;                    // vector of trait b
  
  // pMI (prob mosquito infection)
  int pMI_climate_N;                                // number of observations
  vector[pMI_climate_N] pMI_climate_temp;           // vector of temperatures
  vector[pMI_climate_N] pMI_climate;                // vector of trait pMI

  
  // new climate data for predictions 
  int climate_N_new;                                // number of new temperatures
  vector[climate_N_new] climate_temp_new;           // vector of new temperatures
}

parameters {

  
  // b (prob mosquito infectiousness)
  real<lower=0, upper=0.01> b_climate_constant;     // parameter c
  real<lower=0, upper=24> b_climate_Tmin;           // parameter Tmin
  real<lower=24, upper=45> b_climate_Tmax;          // parameter Tmax
  real<lower=0> b_climate_sigma;                    // noise
  
  // pMI (prob mosquito infection)
  real<lower=0, upper=10> pMI_climate_rmax;          // parameter rmax
  real<lower=10, upper=40> pMI_climate_Topt;        // parameter Topt
  real<lower=0, upper=10> pMI_climate_a;            // parameter a
  real<lower=0> pMI_climate_sigma;                  // noise

}

model {                                             // Fit models to observed data

  
  // b (prob mosquito infectiousness)
  b_climate_constant ~ normal(8.49E-04,0.01);       // prior for c
  b_climate_Tmin ~ normal(17.05,20);                // prior for Tmin
  b_climate_Tmax ~ normal(35.83,20);                // prior for Tmax
  b_climate_sigma ~ uniform(0,100);                 // prior for sigma

  for(j in 1:b_climate_N){                          // gaussian model
    real b_climate_mu = b_climate_constant * b_climate_temp[j] * (b_climate_temp[j] - b_climate_Tmin) * sqrt(b_climate_Tmax - b_climate_temp[j]);
    b_climate[j] ~ normal(b_climate_mu, b_climate_sigma);
  }
  if (b_climate_Tmin > b_climate_Tmax) {
    target += positive_infinity();
  }
  
  // pMI (prob mosquito infection)
  pMI_climate_rmax ~ normal(0.24, 0.03);            // prior for rmax
  pMI_climate_Topt ~ normal(30.08, 0.38);           // prior for Topt
  pMI_climate_a ~ normal(3.60, 0.41);               // prior for a
  pMI_climate_sigma ~ uniform(0,100);               // prior for sigma

  for(k in 1:pMI_climate_N){                        // gaussian model
    real pMI_climate_mu = pMI_climate_rmax * exp(-0.5 * (fabs(pMI_climate_temp[k] - pMI_climate_Topt)/pMI_climate_a)^2);
    pMI_climate[k] ~ normal(pMI_climate_mu, pMI_climate_sigma);
  }

}

generated quantities {
  
  // posterior predictive check values
  real b_climate_ppc[b_climate_N];
  real pMI_climate_ppc[pMI_climate_N];              // pMI (prob mosquito infection)

  vector[climate_N_new] b_climate_new;              // b (prob mosquito infectiousness)
  vector[climate_N_new] pMI_climate_new;            // pMI (prob mosquito infection)

  // ppc b (prob mosquito infectiousness)
  for (q in 1:b_climate_N){
    real b_climate_mu_ppc = b_climate_constant * b_climate_temp[q] * (b_climate_temp[q] - b_climate_Tmin) * sqrt(b_climate_Tmax - b_climate_temp[q]);
    b_climate_ppc[q] = normal_rng(b_climate_mu_ppc, b_climate_sigma);
  }

  // ppc pMI (prob mosquito infection)
  for (r in 1:pMI_climate_N){
    real pMI_climate_mu_ppc = pMI_climate_rmax * exp(-0.5 * (fabs(pMI_climate_temp[r] - pMI_climate_Topt)/pMI_climate_a)^2);;
    pMI_climate_ppc[r] = normal_rng(pMI_climate_mu_ppc, pMI_climate_sigma);
  }


  // predict across new temperatures
  for (xx in 1:climate_N_new){                      
    

    // b (prob mosquito infectiousness)
    if(b_climate_Tmin < climate_temp_new[xx] && b_climate_Tmax > climate_temp_new[xx]){
      b_climate_new[xx] = normal_rng((b_climate_constant * climate_temp_new[xx] * (climate_temp_new[xx] - b_climate_Tmin) * sqrt(b_climate_Tmax - climate_temp_new[xx])), b_climate_sigma);
    }
    else {
      b_climate_new[xx] = 0;
    }

    // pMI (prob mosquito infection)
    pMI_climate_new[xx] = normal_rng(pMI_climate_rmax * exp(-0.5 * (fabs(climate_temp_new[xx] - pMI_climate_Topt)/pMI_climate_a)^2), pMI_climate_sigma);



  }
}