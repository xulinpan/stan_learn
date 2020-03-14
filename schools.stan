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
  int<lower=0> n; //number of schools
  real y[n]; // effect of coaching
  real<lower=0> sigma[n]; // standard errors of effects
}


// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real mu;  // the overall mean effect
  real<lower=0> tau; // the inverse variance of the effect
  vector[n] eta; // standardized school-level effects (see below)
}
transformed parameters {
  vector[n] theta; 
  theta = mu + tau * eta; // find theta from mu, tau, and eta
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  target += normal_lpdf(eta | 0, 1); // eta follows standard normal
  target += normal_lpdf(y | theta, sigma);  // y follows normal with mean theta and sd sigma
}
