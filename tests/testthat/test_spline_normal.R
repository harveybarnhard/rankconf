# library(rstan)
# library(splines)
#
# set.seed(1234)
# num_knots = 10
# spline_degree = 3
# num_basis = num_knots + spline_degree - 1
# X = seq(from=-10, to=10, by=.2)
# knots = unname(quantile(X,probs=seq(from=0, to =1, length.out=num_knots)))
# num_data = length(X)
# a0 = 0.2
# a = rnorm(num_basis, 0, 1)
# B_true = t(bs(X, df=num_basis, degree=spline_degree, intercept = TRUE))
# Y_true = as.vector(a0*X + a%*%B_true)
# Y = Y_true = rnorm(length(X), 0, 0.2)
#
# num_knots = 5
# num_basis = num_knots + spline_degree - 1
# knots = unname(quantile(Y, probs=seq(from=0, to=1, length.out=num_knots)))
# B_mat = t(bs(Y, df=num_basis, degree=spline_degree, intercept=TRUE))
# rstan(options(auto_write = TRUE))
# options(mc.cores = parallel::detectCores())
# spline_model = stan_model(file.path(.github, "rankconf/inst/stan/spline_normal2.stan"))
#
# fit_spline = sampling(
#   spline_model, iter=500, control=list(adapt_delta=0.7),
#   data = list(
#     n = length(Y),
#     y = Y,
#     sigma = 0.9*sqrt(Y^2),
#     num_knots = num_knots,
#     knots = knots,
#     spline_degree = spline_degree,
#     alpha_gamma1 = 0.1,
#     alpha_gamma2 = 0.1
#   ),
#   chains=1,
#   init = list(
#     list(theta = rep(0, 101))
#   ),
#   verbose=TRUE
# )
