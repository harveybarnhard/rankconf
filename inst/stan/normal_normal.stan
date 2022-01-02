data {
  int<lower=0> n;          // number of groups
  real y[n];               // effect estimates
  real<lower=0> sigma[n];  // s.e. of effect estimates
}
parameters {
  real mu;
  real<lower=0> tau;
  vector[n] eta;
}
transformed parameters {
  vector[n] theta;
  theta = mu + tau * eta;
}
model {
  target += normal_lpdf(eta | 0, 1);
  target += normal_lpdf(y | theta, sigma);
}
