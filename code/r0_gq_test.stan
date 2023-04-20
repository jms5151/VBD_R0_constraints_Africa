//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
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
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
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



generated quantities{
  vector[N_new] NmNh_new;                                           // ratio of mosquitoes to humans
  vector[N_new] delta_new;                                          // intrinsic incubation period
  vector[N_new] mu_h_new;                                           // human mortality rate
  vector[N_new] gamma_new;                                          // Human infectivity period
  vector[N_new] alpha_new;                                          // alpha (biting rate)
  vector[N_new] b_new;                                              // b (prob mosquito infectiousness)
  vector[N_new] pMI_new;                                            // pMI (prob mosquito infection)
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

}