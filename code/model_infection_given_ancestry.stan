data {
  int<lower=0> pMI_ancestry_N;                           // number of observations
  int<lower=0> Trials[pMI_ancestry_N];                   // number of trials for each observation
  int<lower=0> Infected[pMI_ancestry_N];                 // number of successes for each observation
  matrix[pMI_ancestry_N, 2] X;                           // matrix of covariates (continuous)
  int<lower=0, upper=1> Virus[pMI_ancestry_N];           // categorical variable
  
  // new data for predictions 
  int pMI_ancestry_N_new;                                // number of new values
  matrix[pMI_ancestry_N_new, 2] X_new;                   // matrix of new of continuous covariates
  int<lower=0, upper=1> Virus_new[pMI_ancestry_N_new];       // new categorical variable

}

parameters {
  real pMI_ancestry_b0;                                  // intercept
  vector[2] beta_ancestry;                               // coefficients for continuous covariates
  // real anc_beta;
  // real logdose_beta;
  real gamma_ancestry;                                   // coefficient for categorical covariate
}

model {
  pMI_ancestry_b0 ~ normal(-13,10);                      // prior for intercept
  beta_ancestry ~ normal(0, 5);                          // priors for continuous variable coefficients
  // anc_beta ~ normal(2.65407, 0.23774);
  // logdose_beta ~ normal(0.78197, 0.04504);
  gamma_ancestry ~ normal(0, 5);                       // prior for categorical variable coefficients
  for (i in 1:pMI_ancestry_N) {
    Infected[i] ~ binomial(Trials[i], inv_logit(pMI_ancestry_b0 + X[i] * beta_ancestry + gamma_ancestry * Virus[i]));
    // Infected[i] ~ binomial(Trials[i], inv_logit(anc_beta * anc[i] + logdose_beta * logdose[i] + gamma_ancestry * Virus[i]));
  }
}


generated quantities {
    // posterior predictive check values
  real<lower=0> Infected_ppc[pMI_ancestry_N];            //
  real<lower=0, upper=1> propInf[pMI_ancestry_N];

  vector[pMI_ancestry_N_new] pMI_ancestry_pred; 

  // ppc 
  for (p in 1:pMI_ancestry_N){
    propInf[p] = inv_logit(pMI_ancestry_b0 + X[p] * beta_ancestry + gamma_ancestry * Virus[p]);
    Infected_ppc[p] = binomial_rng(Trials[p], propInf[p]);
    }
  
  // new predictions
  for(pp in 1:pMI_ancestry_N_new){
    pMI_ancestry_pred[pp] = inv_logit(pMI_ancestry_b0 + X_new[pp] * beta_ancestry + gamma_ancestry * Virus_new[pp]);
  }
  
}
