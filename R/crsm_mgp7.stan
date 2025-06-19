// adapted from crsm_mgp6.stan, updated so that the prior for mu can be set by user rather than hard-coded.
functions {
  // competing risks cox proportional hazards model, log likelihood
  real crcox_lpdf(vector y, matrix Xbeta, matrix delta, matrix lambda_cal, matrix H) {
    return sum(delta .* (log(lambda_cal) + Xbeta) - exp(Xbeta) .* H);
  }
}

data {
  int<lower=1> n; // number of data points
  int<lower=1> p; // number of covariates
  int<lower=1> k; // number of bins for piecewise linear baseline hazard rates
  int<lower=1> m; // number of risk types
  vector[n] y; // time to event data
  vector[k+1] s; // knots for the piecewise linear baseline hazard rates
  matrix[n,m] delta; // indicator of risk type
  matrix[n,p] X; // design matrix
  real<lower=0> kappa0; // for prior of psi0
  real<lower=0> a; // for prior of kappa1
  real<lower=0> b; // for prior of kappa1
  real mu0; // for prior of mu
  real<lower=0> tau0; // for prior of mu
}

transformed data {
  array[n] int<lower=1,upper=k> ki; // bin for each event
  matrix[k,m] s_diff = rep_matrix(s[2:(k+1)]-s[1:k],m); // bin widths, matrix form for calculation
  matrix[n,m] T; // time since last knot, matrix form for calculation
  row_vector[m] zeros = zeros_row_vector(m);

  for (i in 1:n) {
    for (l in 1:k){
      if (y[i]>s[l] && y[i] <= s[l+1]){
        ki[i] = l;
        T[i] = rep_row_vector(y[i] - s[l],m);
      }
    }
  }
}

parameters {
  matrix[p,m] beta;
  row_vector<lower=0>[m] psi0;
  matrix<lower=0>[k-1,m] psi;
  vector[p] mu;
  real<lower=0> s2;
  real<lower=0> kappa1; // for multiplicative gamma process
}

transformed parameters {
  matrix[k,m] lambda; // baseline hazard rates
  
  lambda[1] = psi0;
  for (l in 2:k){
    lambda[l] = psi[l-1] .* lambda[l-1];
  }
}

model {
  matrix[n,m] H; // baseline cumulative hazard rates for each ti
  matrix[n,m] lambda_cal; // lambda for each ti
  matrix[k,m] cum_lambda1;
  matrix[k,m] cum_lambda2;
  matrix[p,m] mu_rep;
  matrix[n,m] Xbeta;
  vector[(k-1)*m] psi_1d;
  vector[p*m] mu_1d;
  vector[p*m] beta_1d;

  // get lambda_cal
  lambda_cal = lambda[ki];
  
  // get H
  cum_lambda1 = lambda .* s_diff;
  cum_lambda2[1] = zeros;
  for (l in 2:k){
    cum_lambda2[l] = cum_lambda2[l-1] + cum_lambda1[l-1];
  }
  H = cum_lambda2[ki];

  H = H + lambda_cal .* T;
  mu_rep = rep_matrix(mu,m);
  mu_1d = to_vector(mu_rep);
  beta_1d = to_vector(beta);
  psi_1d = to_vector(psi);
  Xbeta = X * beta;

  target += crcox_lupdf(y | Xbeta, delta, lambda_cal, H);
  target += inv_gamma_lupdf(s2 |1,1);
  target += normal_lupdf(beta_1d | mu_1d, sqrt(s2));
  target += normal_lupdf(mu | mu0, tau0);
  target += gamma_lupdf(psi0 | kappa0, kappa0);
  target += gamma_lupdf(psi_1d | kappa1, kappa1);
  target += inv_gamma_lupdf(kappa1 | a, b);
}
