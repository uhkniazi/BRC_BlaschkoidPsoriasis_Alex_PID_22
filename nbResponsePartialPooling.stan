data {
  int<lower=1> Ntotal; // number of observations
  int<lower=1> Ncol; // total number of columns in model matrix
  int<lower=1> NscaleBatches;
  matrix[Ntotal, Ncol] X; // model matrix
  int y[Ntotal]; // response variable nb distributed
  int<lower=1, upper=NscaleBatches> NBatchMap[Ncol];
}
// transformed data {
  // }
parameters {
  // parameters to estimate in the model
  vector[Ncol] betas; // regression parameters
  real populationMean;//[NscaleBatches];
  real<lower=0.01> sigmaRan[NscaleBatches]; // group level errors
  real<lower=0> phi_scaled; // over dispersion on square root scale 
}
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  real phi;
  phi = phi_scaled^2;
  // fitted value
  mu = X * betas;
  mu = exp(mu + populationMean);
}
model {
  real sigmaRan_expanded[Ncol];
  // using diffuse prior
  sigmaRan ~ exponential(1);
  populationMean ~ cauchy(0, 1);
  phi_scaled ~ normal(0, 10); // weak prior on square root scale
  // vector expansion by mapping to a larger vector/array
  sigmaRan_expanded = sigmaRan[NBatchMap];
  //populationMean_expanded = populationMean[NBatchMap];
  betas ~ normal(0, sigmaRan_expanded); //prior for the betas
  // likelihood function
  y ~ neg_binomial_2(mu, phi);
}
