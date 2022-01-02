// Modified from Stan documentation, Milad Kharratzadeh:
// https://mc-stan.org/users/documentation/case-studies/splines_in_stan.html
functions {
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order);
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order) {
    // INPUTS:
    //    t:          the points at which the b_spline is calculated
    //    ext_knots:  the set of extended knots
    //    ind:        the index of the b_spline
    //    order:      the order of the b-spline
    vector[size(t)] b_spline;
    vector[size(t)] w1 = rep_vector(0, size(t));
    vector[size(t)] w2 = rep_vector(0, size(t));
    if(order==1)
      for (i in 1:size(t)) // B-splines of order 1 are piece-wise constant
        b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind+1]);
    else {
      if (ext_knots[ind] != ext_knots[ind+order-1])
        w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t))) /
             (ext_knots[ind+order-1] - ext_knots[ind]);
      if (ext_knots[ind+1] != ext_knots[ind+order])
        w2 = 1 - (to_vector(t) - rep_vector(ext_knots[ind+1], size(t))) /
                 (ext_knots[ind+order] - ext_knots[ind+1]);
      // Calculating the B-spline recursively as linear interpolation of two
      // lower-order splines
      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order-1) +
                 w2 .* build_b_spline(t, ext_knots, ind+1, order-1);
    }
    return b_spline;
  }
}


data {
  int<lower=0> n;          // number of estimates
  real y[n];               // point estimates
  real<lower=0> sigma[n];  // s.e. of point estimates
  int<lower=0> num_knots;  // number of knots
  vector[num_knots] knots; // sequence of knots
  int spline_degree;       // the degree of spline (equal to order - 1)
  real alpha_gamma1;       // First inverse gamma hyperprior for spline
  real alpha_gamma2;       // Second inverse gamma hyperprior for spline
}

transformed data {
  // Declare variables
  int num_basis = num_knots + spline_degree - 1; // total number of B-splines
  matrix[num_basis, num_basis] P; // Second difference matrix
  vector[spline_degree + num_knots] ext_knots_temp;
  vector[2*spline_degree + num_knots] ext_knots; // set of extended knots

  // Second difference matrix
  P = rep_matrix(0, num_basis, num_basis);
  for (i in 1:num_basis) {
    P[i,i] = 2;
    if(i > 1)
      P[i,i-1] = -1;
    if(i < num_basis)
      P[i, i+1] = -1;
  }

  // Knots for density function
  ext_knots_temp = append_row(rep_vector(knots[1], spline_degree), knots);
  ext_knots = append_row(ext_knots_temp, rep_vector(knots[num_knots], spline_degree));
}

parameters {
  real<lower=0> sigma_gamma;
  vector[num_basis] gamma;
  real theta[n];
}

transformed parameters {
  matrix[n, num_basis] B;
  vector[num_basis] expgamma;

  // Matrix of B-splines (evaluation points change at each iteration)
  for (i in 1:num_basis) {
    B[:,i] = build_b_spline(theta, to_array_1d(ext_knots), i, spline_degree + 1);
  }
  B[n, num_basis] = 1;

  // Exponentiated coefficients
  expgamma = exp(gamma);
}

model {
  // Priors on distribution of spline coefficients
  target += inv_gamma_lpdf(sigma_gamma | alpha_gamma1, alpha_gamma2);
  target += -num_basis*0.5*log(sigma_gamma) - quad_form(P, gamma)*0.5/sigma_gamma;

  // Likelihood of theta given spline distribution
  target += sum(log(B*expgamma)) - n*log(sum(expgamma));

  // Likelihood of observations given theta
  target += normal_lpdf(y | theta, sigma);
}

generated quantities {
  vector[n] density;
  density = B*expgamma / sum(expgamma);
}
