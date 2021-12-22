#' Construct confidence intervals around ranks.
#'
#' @param y estimates
#' @param sig2 variance of estimates
#' @param type What error should be controlled?
#' @param method What method should be used to control error rates?
#' @param alpha The rate at which the error should be controlled. 0.05 by default.
#' @param k Number of errors to control, used for familywise error rates. 1 by default.
#' @param best "min" for rank 1 to correspond to lowest estimate.
#' "max" for rank 1 to correspond to greatest estimate.
#' @param thr Number of threads to use for parallel processing. 2 by default.
#' @param nchains For Bayesian Analysis: How many markov chains?
#' @param nwarmup For Bayesian Analysis: How many warmup iterations?
#' @param niter   For Bayesian Analysis: How many total iterations?
#' @export
rankconf = function(y,
                    sig2,
                    type="PCER",
                    method="NAIVE",
                    alpha=0.05,
                    k=1,
                    best="min",
                    thr=2,
                    nchains=4,
                    nwarmup=1500,
                    niter=2500
                    ){
  # Handle Input ===============================================================
  if(any(is.na(y) | is.na(sig2) | sig2 <= 0)){
    stop("y and sig2 must not contain NA values, and sig2 must be positive")
  }

  # Make sure error types exist and an appropriate method is called
  type   = toupper(type)
  method = toupper(method)
  if(type=="FDR"){
    if(!method%in%c("BH", "BY")){
      stop("For FDR, the following methods are allowed: Benjamini-Hochberg (BH) or Benjamini-Yekutieli (BY)")
    }
  }
  else if(type=="FWER"){
    if(!method%in%c("B", "HB", "R")){
      stop("For FWER, the following methods are allowed: Bonferroni (B), Holm-Bonferroni (HB), or Romano re-sampling (R)")
    }
  }

  # Check to make sure inputs make sense
  n = length(y)
  if(length(sig2)!=n){
    stop("y and sig2 must be of same length")
  }
  if(best=="min"){
    message("Rank 1 corresponds to the lowest estimate\n")
  }else if(best=="max"){
    message("Rank 1 corresponds to the highest estimate\n")
    y = -y
  }

  # Frequentist multiple testing methods
  if(type%in%c("PCER","FDR", "FWER")) {
    output = rankconf_multitest(n, y, sig2, type, method, alpha, k, thr)
  }
  # Bayesian posterior inference methods
  else if(type%in%c("BAYES")) {
    output = rankconf_bayes(n, y, sig2, type, method, alpha, thr, nchains, nwarmup, niter)
  }

  # Add more data to output
  output[["data"]]$y = ifelse(best=="max", -y, y)
  output[["data"]]$sig2 = sig2
  output[["data"]]$y_rank = data.table::frank(y, ties.method="min", na.last=T)
  return(output)
}

# Frequentist multiple testing methods =========================================
rankconf_multitest = function(n, y, sig2, type, method, alpha, k, thr) {
  # Calculate naive one-sided p-values of all differences
  diffmat = matrix(selfouter(y, '-')/sqrt(selfouter(sig2, "+")), n, n)
  diag(diffmat) = NA
  if(!(type=="FWER" & method=="R")){
    pvals = 1-stats::pnorm(diffmat)
  }else{
    pvals = NA
    diffmat = abs(diffmat)
  }
  gc()

  # Find which tests to reject using the given method
  reject = do.call(
    paste(c("rej", type, method), collapse="_"),
    list(
      diffmat=diffmat,
      pvals=pvals,
      alpha=alpha,
      k=k,
      sig2=sig2,
      y=y,
      thr=thr
    )
  )
  diag(reject) = FALSE

  # Create a "reject and positive difference" matrix.
  # A TRUE in the [i,j]th index indicates that the ith observation is
  # significantly better than the jth observation.
  if(type=="FWER" & method=="R"){
    diffmat = matrix(selfouter(y, '-')/sqrt(selfouter(sig2, "+")), n, n)
  }
  diffmat = (diffmat > 0) & reject & !is.na(diffmat)
  rm(reject)
  gc()
  # Calculate number of ranks above and below.
  # Summing over the jth column is equivalent to adding up the number of
  # observations that are significantly worse than the jth observation.
  # Summing over the ith row is equivalent to adding up the number of
  # observations that are significantly better than the ith observation.
  out_data = data.table::data.table(
    L = rowSums(diffmat) + 1,
    U = n - colSums(diffmat)
  )
  return(list(data = out_data))
}

# Bayesian posterior inference methods =========================================
rankconf_bayes = function(n, y, sig2, type, method, alpha, thr, nchains, nwarmup, niter) {
  # Format the data for Stan
  stan_data = list(
    J = n,
    y = y,
    sigma = sqrt(sig2)
  )

  # Sample from posterior draws
  stan_fit <- rstan::sampling(
    stanmodels$normal_normal,    # Stan program
    data = stan_data,             # named list of data
    chains = nchains,            # number of Markov chains
    warmup = nwarmup,            # number of warmup iterations per chain
    iter = niter,                # total number of iterations per chain
    cores = thr,                 # number of cores (could use one per chain)
    refresh = 1                  # show progress
  )

  # Extract draws for the group-level means
  ests = rstan::extract(stan_fit)$theta

  # Construct all pairwise differences in each posterior draw
  diffs = matrix(0, nrow=n, ncol=n)
  pb = utils::txtProgressBar(min = 1, max = n, initial = 1, style=3)
  for(i in 1:n){
    utils::setTxtProgressBar(pb,i)
    diffs = ((i-1)*diffs + outer(ests[i,], ests[i,], ">"))/i
  }
  close(pb)

  # Construct the rank credible sets
  diffs = diffs > 1-alpha
  out_data = data.table::data.table(
    L = rowSums(diffs) + 1,
    U = nrow(diffs) - colSums(diffs)
  )
  return(list(
    data   = out_data,
    est    = ests
  ))
}
