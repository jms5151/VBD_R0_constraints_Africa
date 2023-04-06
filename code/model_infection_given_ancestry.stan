data {
  int<lower=0> pMI_ancestry_N;                           // number of observations
  int<lower=0> pMI_ancestry_Trials[pMI_ancestry_N];                   // number of trials for each observation
  int<lower=0> pMI_ancestry_Infected[pMI_ancestry_N];                 // number of successes for each observation
  matrix[pMI_ancestry_N, 2] pMI_ancestry_X;                           // matrix of covariates (continuous)
  int<lower=0, upper=1> pMI_ancestry_Virus[pMI_ancestry_N];           // categorical variable
  
  // new data for predictions 
  int pMI_ancestry_N_new;                                // number of new values
  matrix[pMI_ancestry_N_new, 2] pMI_ancestry_X_new;                   // matrix of new of continuous covariates
  int<lower=0, upper=1> pMI_ancestry_Virus_new[pMI_ancestry_N_new];       // new categorical variable

}

parameters {
  real pMI_ancestry_b0;                                  // intercept
  vector[2] pMI_ancestry_beta;                               // coefficients for continuous covariates
  real pMI_ancestry_gamma;                                   // coefficient for categorical covariate
}

model {
  pMI_ancestry_b0 ~ normal(-13,10);                      // prior for intercept
  pMI_ancestry_beta ~ normal(0, 5);                          // priors for continuous variable coefficients
  pMI_ancestry_gamma ~ normal(0, 5);                       // prior for categorical variable coefficients
  
  for (i in 1:pMI_ancestry_N) {
    pMI_ancestry_Infected[i] ~ binomial(pMI_ancestry_Trials[i], inv_logit(pMI_ancestry_b0 + pMI_ancestry_X[i] * pMI_ancestry_beta + pMI_ancestry_gamma * pMI_ancestry_Virus[i]));
  }
}


generated quantities {
    // posterior predictive check values
  real<lower=0> pMI_ancestry_Infected_ppc[pMI_ancestry_N];            //
  real<lower=0, upper=1> pMI_ancestry_propInf[pMI_ancestry_N];

  vector[pMI_ancestry_N_new] pMI_ancestry_pred; 

  // ppc 
  for (p in 1:pMI_ancestry_N){
    pMI_ancestry_propInf[p] = inv_logit(pMI_ancestry_b0 + X[p] * pMI_ancestry_beta + pMI_ancestry_gamma * pMI_ancestry_Virus[p]);
    pMI_ancestry_Infected_ppc[p] = binomial_rng(pMI_ancestry_Trials[p], pMI_ancestry_propInf[p]);
    }
  
  // new predictions
  for(pp in 1:pMI_ancestry_N_new){
    pMI_ancestry_pred[pp] = inv_logit(pMI_ancestry_b0 + pMI_ancestry_X_new[pp] * pMI_ancestry_beta + pMI_ancestry_gamma * pMI_ancestry_Virus_new[pp]);
  }
  
}
