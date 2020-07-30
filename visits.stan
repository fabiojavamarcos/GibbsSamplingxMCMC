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
  int<lower=1> N;
  int<lower=1> p;
  matrix[N,p] X;
  vector[N] y;
  vector[p] m0;
  matrix[p,p] prec0;
  real<lower=0> nu0;
  real<lower=0> sig20;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
   vector[p] beta;
  real<lower=0> sigma2;
}

transformed parameters{
  real<lower=0> sigma=sqrt(sigma2);
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  sigma2 ~ scaled_inv_chi_square(nu0, sig20);
  y ~ normal(X*beta, sigma);
  beta ~ multi_normal_prec(m0, prec0);
}

