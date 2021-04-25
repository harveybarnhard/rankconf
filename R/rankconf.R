#' Construct confidence intervals around ranks.
#'
#' @param y estimates
#' @param sig2 variance of estimates
#' @param type What error should be controlled?
#' @param alpha The rate at which the error should be controlled. 0.05 by default
#' @param k For k-FWER, number of errors to be bounded
#' @param best "min" for rank 1 to correspond to lowest estimate.
#' "max" for rank 1 to correspond to greatest estimate.
#' @param thr Number of threads to utilize. Number of cores minus one by default
#' @export
rankconf = function(y,
                    sig2,
                    type="FDR",
                    alpha=0.05,
                    k=NA,
                    best="min",
                    thr=parallel::detectCores()-1,
                    verbose=F){
  # Handle Input ===============================================================
  if(any(is.na(y) | is.na(sig2) | sig2 <= 0)){
    stop("y and sig2 must not contain NA values, and sig2 must be positive")
  }

  # Make sure error types exist
  # TODO: Separate this between error rates and methods to control error rates
  type = toupper(type)

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

  # Calculate naive one-sided p-values of all differences
  diffmat = matrix(selfouter(y, '-')/sqrt(selfouter(sig2, "+")), n, n)
  if(type!="BKFWER"){
    diffmat = matrix(selfouter(y, '-')/sqrt(selfouter(sig2, "+")), n, n)
    pvals = 1-pnorm(diffmat)
  }else{
    pvals = NA
    diffmat = abs(matrix(selfouter(y, '-')/sqrt(selfouter(sig2, "+")), n, n))
    diag(diffmat) = NA
  }
  gc()

  # Find which tests to reject using the given method
  reject = do.call(paste0("rej", type), list(diffmat=diffmat,
                                             pvals=pvals,
                                             alpha=alpha,
                                             k=k,
                                             sig2=sig2,
                                             y=y,
                                             thr=thr))
  diag(reject) = FALSE

  # Create a "reject and positive difference" matrix.
  # A TRUE in the [i,j]th index indicates that the ith observation is
  # significantly better than the jth observation.
  if(type=="BKFWER"){
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
  lowerrank = rowSums(diffmat) + 1
  upperrank = n - colSums(diffmat)
  return(list(
    L = lowerrank,
    U = upperrank
  ))
}

# Naive method (PCER)
rejNAIVE = function(pvals, alpha, ...){
  return(pvals < alpha/2)
}

# Bonferroni (FWER)
rejBONF = function(pvals, alpha, ...){
  n = nrow(pvals)
  return(pvals < alpha/(n^2-n))
}

# Holm-Bonferroni (FWER)
rejFWER = function(pvals, alpha, ...){
  n = nrow(pvals)
  m = n^2 - n
  # Rank pvals from lowest to highest
  pranks = matrix(data.table::frank(c(pvals), ties.method="min", na.last=T), n, n)

  # Find first p value not low enough for rejection, reject everything below
  min_rank = min(c(pranks[pvals>=alpha/(m +1-pranks)], Inf), na.rm=T)
  return(pranks < min_rank)
}

# Benjamini-Yekutieli (FDR)
rejFDR = function(pvals, alpha, ...){
  n = nrow(pvals)
  m = n^2 - n
  cm = log(m) - digamma(1) + 1/(2*m)
  pranks = matrix(data.table::frank(c(pvals), ties.method="min", na.last=T), n, n)
  max_rank = max(c(pranks[pvals<=pranks*alpha/(m*cm)], -Inf), na.rm=T)
  return(pranks < max_rank)
}

# Bonferroni (KBONF)
rejKBONF = function(pvals, alpha, k, ...){
  n = nrow(pvals)
  m = n^2 - n
  return(pvals < k*alpha/m)
}

# Holm-Bonferroni (KFWER)
rejKFWER = function(pvals, alpha, k, ...){
  n = nrow(pvals)
  m = n^2 - n
  # Rank pvals from lowest to highest
  pranks = matrix(data.table::frank(c(pvals), ties.method="min", na.last=T), n, n)

  # Create values against which to compare p values
  compvals = matrix(NA, n, n)
  ind = pranks<=k
  compvals[ind]  = k*alpha/m
  compvals[!ind] = k*alpha/(m+k-pranks[!ind])
  diag(compvals) = NA

  # Find first p value not low enough for rejection, reject everything below
  min_rank = min(pranks[pvals>=compvals], na.rm=T)
  return(pranks < min_rank)
}

# (REKFWER) https://arxiv.org/pdf/0710.2258.pdf Algorithm 2.2 ==================

# Function that returns one bootstrap sample of the kth largest value
sampfun = function(sigmat, k, distfun, ...){
  booty = do.call(distfun, list(...))
  return(
    kmax(
      abs(selfouter(booty, "-")/sigmat),
      k
    )
  )
}

# Main function
rejBKFWER = function(diffmat, sig2, alpha, k, R=1000, distfun="rnorm", thr=0, ...){
  # Parse arguments
  if(distfun=="rnorm"){
    distfun = function(sd){
      Rfast::Rnorm(length(sd), m=0, s=1)*sd
    }
  }

  # Initialize cluster for parallel processing
  cl = parallel::makeCluster(thr)
  parallel::clusterExport(
    cl, varlist=c("kmax", "selfouter"),
    envir=environment()
  )
  on.exit(parallel::stopCluster(cl))

  # Initialize rejection matrix with no rejections
  reject = matrix(F, nrow=length(sig2), ncol=length(sig2))

  # Calculate size of sample to include in each resample
  n = nrow(diffmat)
  numind = n^2 - n

  #  Create covariance matrix of test statistics
  sigmat = sqrt(matrix(selfouter(sig2, FUN = "+"), n, n))
  diag(sigmat) = -diag(sigmat)

  # Step 1, calculate 1-alpha critical values
  #     a) Reject all values greater than critical value
  #     b) If k or fewer observations are rejected, then stop
  s = sqrt(sig2)
  kmaxdist = parallel::parSapply(
    cl=cl,
    X=rep(list(1), R),
    FUN=function(x) do.call(
      sampfun, list(sigmat=sigmat, k=k, distfun=distfun, sd=s)
    )
  )
  # Return the quantile of the distribution
  crit = quantile(kmaxdist, probs=1-alpha, names=FALSE)
  numrej = rejupdate(reject, diffmat, crit)
  cat("\r Step 1. ", numrej, " new rejections.")

  # Step 2,3,... only if there were at least k rejections in the first step
  newrej = numrej >= k
  j = 2
  while(newrej){
    # Calculate critical value over unrejected tests and k-1 least significant
    # tests that have already been rejected
    if(k > 1) {
      c = kmin(x=diffmat[reject], k=k-1)
    }else {
      c = NA
    }
    numind = sigupdate(reject, sigmat, diffmat, c, numind, k)
    kmaxdist = parallel::parSapply(
      cl=cl,
      X=rep(list(1), R),
      FUN=function(x) do.call(
        sampfun, list(sigmat=sigmat, k=k, distfun=distfun, sd=s)
      )
    )

    # Find the 1-alpha quantile of the distribution
    crit = quantile(kmaxdist, probs=1-alpha, names=FALSE)

    # Update the rejection matrix based on the 1-alpha critical value
    # Continue if there are any new rejections
    numrej = rejupdate(reject, diffmat, crit)
    newrej = numrej > 0

    # Display status
    cat("\r Step", j,". ", numrej, " new rejections.")
    flush.console()
    j = j + 1
    gc()
  }
  return(reject)
}
