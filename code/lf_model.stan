data {
  
  // lifespan (1/mosquito mortality rate)
  int lf_climate_N;                                             // number of observations
  vector[lf_climate_N] lf_climate_temp;                         // vector of temperatures
  vector[lf_climate_N] lf_climate;                              // vector of trait lf

  // new data
  int N_new;
  vector[N_new] temp_new;
}

parameters {
  // lifespan (1/mosquito mortality rate)
  real<lower=-0.5, upper=0.5> lf_climate_constant;              // parameter c
  real<lower=0, upper=24> lf_climate_Tmin;                      // parameter Tmin
  real<lower=24, upper=45> lf_climate_Tmax;                     // parameter Tmax
  real<lower=0> lf_climate_sigma;                               // noise

  // guassian
  // real<lower=0, upper=10> lf_climate_rmax;                     // parameter rmax
  // real<lower=10, upper=40> lf_climate_Topt;                    // parameter Topt
  // real<lower=0, upper=10> lf_climate_a;                        // parameter a
  // real<lower=0> lf_climate_sigma;                              // noise


}

model {                                                         // Fit models to observed data
  // lifespan (1/mosquito mortality rate)
  lf_climate_constant ~ normal(-1.48E-01,1);                    // prior for c
  lf_climate_Tmin ~ normal(9.16,1);                             // prior for Tmin
  lf_climate_Tmax ~ normal(37.73,1);                            // prior for Tmax
  lf_climate_sigma ~ normal(0.01,0.1);                          // prior for sigma

  // quadratic
  // for(m in 1:lf_climate_N){
  //   real lf_climate_mu = lf_climate_constant * (lf_climate_temp[m] - lf_climate_Tmax) * (lf_climate_temp[m] - lf_climate_Tmin);
  //   lf_climate[m] ~ normal(lf_climate_mu, lf_climate_sigma);
  // }
  // if (lf_climate_Tmin > lf_climate_Tmax) {
  //   target += positive_infinity();
  // }
  
  // briere
    for(m in 1:lf_climate_N){
      real lf_climate_mu = lf_climate_constant * lf_climate_temp[m] * (lf_climate_temp[m] - lf_climate_Tmin) * sqrt(lf_climate_Tmax - lf_climate_temp[m]);
      lf_climate[m] ~ normal(lf_climate_mu, lf_climate_sigma);
  }
  if (lf_climate_Tmin > lf_climate_Tmax) {
    target += positive_infinity();
  }

  // // guassian
  // lf_climate_rmax ~ normal(0.24, 1);                            // prior for rmax
  // lf_climate_Topt ~ normal(30.08, 1);                           // prior for Topt
  // lf_climate_a ~ normal(3.60, 1);                               // prior for a
  // lf_climate_sigma ~ normal(0.01,0.1);                          // prior for sigma
  // 
  // for(k in 1:lf_climate_N){                        
  //   real lf_climate_mu = lf_climate_rmax * exp(-0.5 * (fabs(lf_climate_temp[k] - lf_climate_Topt)/lf_climate_a)^2);
  //   lf_climate[k] ~ normal(lf_climate_mu, lf_climate_sigma);
  // }
  
}


generated quantities {
  real lf_climate_ppc[lf_climate_N];                              // lifespan (1/mosquito mortality rate)
  vector[N_new] lf_climate_new;                                   // lifespan (1/mosquito mortality rate)

  // ppc lifespan (1/mosquito mortality rate)
  // quadratic
  // for (u in 1:lf_climate_N){
  //   real lf_climate_mu_ppc = lf_climate_constant * (lf_climate_temp[u] - lf_climate_Tmax) * (lf_climate_temp[u] - lf_climate_Tmin);
  //   lf_climate_ppc[u] = normal_rng(lf_climate_mu_ppc, lf_climate_sigma);
  // }
  
  // briere
  for (p in 1:lf_climate_N){
    real lf_climate_mu_ppc = lf_climate_constant * lf_climate_temp[p] * (lf_climate_temp[p] - lf_climate_Tmin) * sqrt(lf_climate_Tmax - lf_climate_temp[p]);
    lf_climate_ppc[p] = normal_rng(lf_climate_mu_ppc, lf_climate_sigma);
  }

  // ppc guassian
  // for (r in 1:lf_climate_N){
  //   real lf_climate_mu_ppc = lf_climate_rmax * exp(-0.5 * (fabs(lf_climate_temp[r] - lf_climate_Topt)/lf_climate_a)^2);
  //   lf_climate_ppc[r] = normal_rng(lf_climate_mu_ppc, lf_climate_sigma);
  // }

  // predict on new data
  for(zz in 1:N_new){
    
    // lifespan (1/mosquito mortality rate)
    // quadratic
    // if(lf_climate_Tmin < temp_new[zz] && lf_climate_Tmax > temp_new[zz]){
    //   lf_climate_new[zz] = normal_rng((lf_climate_constant * (temp_new[zz] - lf_climate_Tmax) * (temp_new[zz] - lf_climate_Tmin)), lf_climate_sigma);
    // }
    // else {
    //   lf_climate_new[zz] = 0;
    
    // briere
    if(lf_climate_Tmin < temp_new[zz] && lf_climate_Tmax > temp_new[zz]){
      lf_climate_new[zz] = normal_rng((lf_climate_constant * temp_new[zz] * (temp_new[zz] - lf_climate_Tmin) * sqrt(lf_climate_Tmax - temp_new[zz])), lf_climate_sigma);
    }
    else {
      lf_climate_new[zz] = 0;
    }

    // guassian
    // lf_climate_new[zz] = normal_rng(lf_climate_rmax * exp(-0.5 * (fabs(temp_new[zz] - lf_climate_Topt)/lf_climate_a)^2), lf_climate_sigma);

    }
    
}

