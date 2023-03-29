// not working, nor is example model: https://mc-stan.org/docs/2_21/stan-users-guide/logistic-probit-regression-section.html

data {
 int<lower=0> pMI_ancestry_N;                            // number of observations
 vector[pMI_ancestry_N] pMI_ancestry;               // vector of ancestry levels
//  // vector[pMI_ancestry_N] virus;                  // vector of virus type
//  // vector[pMI_ancestry_N] dose;                   // vector of dose levels
 int<lower=0, upper=1> pMI_ancestry_pred[pMI_ancestry_N];           // integer vector of pred pMI | ancestry
// 
//  // int ancestry_N_new;                            // number of new observations to predict on
//  // vector[ancestry_N_new] aa_ancestry_new;        // vector of proportion aa ancestry 
}

parameters {
  real pMI_ancestry_intercept;
  real beta; // vector of coeffients
}

model {
  pMI_ancestry_pred ~ bernoulli_logit(pMI_ancestry_intercept + beta * pMI_ancestry);
}

// generated quantities {
//   //ppc
//   int pMI_ancestry_ppc[pMI_ancestry_N];
// 
//   pMI_ancestry_ppc ~ bernoulli_logit_rng(pMI_ancestry_intercept + beta * ancestry);
// }