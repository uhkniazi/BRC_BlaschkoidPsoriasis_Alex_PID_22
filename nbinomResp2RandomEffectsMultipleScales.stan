data {
  int<lower=1> Ntotal; // number of observations
  int<lower=1> Nclusters1; // number of levels for group 1 for random intercepts
  int<lower=1> Nclusters2; // number of levels for group 2 for random intercepts
  int<lower=1> NScaleBatches1; // number of batches of scale terms for group 1 (subgroups in group 1) 
  int<lower=1> NScaleBatches2; // number of batches of scale terms for group 2 (subgroups in group 2) 
  int<lower=1, upper=Nclusters1> NgroupMap1[Ntotal]; // mapping variable to map each group 1 parameter to a data point 
  int<lower=1, upper=Nclusters2> NgroupMap2[Ntotal]; // mapping variable to map each group 2 parameter to a data point 
  int<lower=1, upper=NScaleBatches1> NBatchMap1[Nclusters1]; // expanding vector to map each scale term to the relevant 
                                                    // set of coefficients from group 1
  
  int<lower=1, upper=NScaleBatches2> NBatchMap2[Nclusters2]; // expanding vector to map each scale term to the relevant 
                                                    // set of coefficients from group 2
  
  int y[Ntotal]; // response variable
  // additional parameters
  //real gammaShape; // hyperparameters for the gamma distribution for batches of scales 
  //real gammaRate;
  real intercept;
  real intercept_sd;
  int<lower=1> Nphi; // number of phi terms for each subset of observations 
  int<lower=1, upper=Nphi> NphiMap[Ntotal]; // mapping variable 
}
// transformed data {
  // }
parameters {
  // parameters to estimate in the model
  real betas; // constant intercept term
  real<lower=0.01> sigmaRan1[NScaleBatches1]; // random effect standard deviations for sub-batches in group 1
  real<lower=0.01> sigmaRan2[NScaleBatches2]; // random effect standard deviations for sub-batches in group 2
  vector[Nclusters1] rGroupsJitter1; // number of random jitters for each level of cluster/group 1
  vector[Nclusters2] rGroupsJitter2; // number of random jitters for each level of cluster/group 2
  real<lower=0> phi_scaled[Nphi]; // over dispersion on square root scale 
}
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  real phi[Nphi];
  for (i in 1:Nphi) {
    phi[i] = phi_scaled[i]^2;
    }
  mu = exp(betas + rGroupsJitter1[NgroupMap1] + rGroupsJitter2[NgroupMap2]);
}
model {
  real sigmaRan1_expanded[Nclusters1]; 
  real sigmaRan2_expanded[Nclusters2];
  //sigmaRan1 ~ gamma(gammaShape, gammaRate);
  sigmaRan1 ~ exponential(1);
  sigmaRan2 ~ exponential(1);
  betas ~ normal(intercept, intercept_sd);
  phi_scaled ~ normal(0, 10); // weak prior on square root scale
  
  // random effects sample
  sigmaRan1_expanded = sigmaRan1[NBatchMap1];
  rGroupsJitter1 ~ normal(0, sigmaRan1_expanded);
  
  sigmaRan2_expanded = sigmaRan2[NBatchMap2];
  rGroupsJitter2 ~ normal(0, sigmaRan2_expanded);
  
  // likelihood function
  y ~ neg_binomial_2(mu, phi[NphiMap]);
}
