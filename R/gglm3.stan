data {
  int<lower=1> n;
  int<lower=1> p;
  vector[n] y;
  matrix[n, p] X;
  real mu0;
  real<lower=0> tau2;
  real<lower=0> a1;
  real<lower=0> b1;
  real<lower=0> a2;
  real<lower=0> b2;
}

parameters {
  vector[p] beta;
  real<lower=0> s2;
  real<lower=0> alpha;
  real mu;
}

transformed parameters {
  vector<lower=0>[n] m;
  m = exp(X * beta);
}

model {
  vector[n] alphas = rep_vector(alpha,n);
  target += gamma_lpdf(y | alpha, alphas./m);
  target += normal_lpdf(beta | mu, sqrt(s2));
  target += normal_lpdf(mu | mu0, sqrt(tau2));
  if (a1==0) {
    target += cauchy_lpdf(s2 |0,1);
  } else {
    target += inv_gamma_lpdf(s2 | a1, b1);
  }
  if (a2==0) {
    target += cauchy_lpdf(alpha |0,1);
  } else {
    target += gamma_lpdf(alpha | a2, b2);
  }
}

