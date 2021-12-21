# Bonferroni (FWER) ============================================================
rej_FWER_B = function(pvals, alpha, k, ...){
  n = nrow(pvals)
  return(pvals < k*alpha/(n^2-n))
}

# Holm-Bonferroni (FWER) =======================================================
rej_FWER_HB = function(pvals, alpha, k, ...){
  n = nrow(pvals)
  m = n^2 - n
  # Rank pvals from lowest to highest
  pranks = matrix(data.table::frank(c(pvals), ties.method="min", na.last=T), n, n)

  # Find first p value not low enough for rejection, reject everything below
  if(k==1){
    compvals = m + 1 - pranks
  }else if(k>1){
    # Create values against which to compare p values
    compvals = matrix(NA, n, n)
    ind = pranks<=k
    compvals[ind]  = k*alpha/m
    compvals[!ind] = k*alpha/(m+k-pranks[!ind])
    diag(compvals) = NA
  }
  min_rank = min(c(pranks[pvals>=compvals], Inf), na.rm=T)
  return(pranks < min_rank)
}

# Romano re-sampling, https://arxiv.org/pdf/0710.2258.pdf Algorithm 2.2 ========

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
rej_FWER_R = function(diffmat, sig2, alpha, k, R=1000, distfun="rnorm", thr=0, ...){
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
  crit = stats::quantile(kmaxdist, probs=1-alpha, names=FALSE)
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
    crit = stats::quantile(kmaxdist, probs=1-alpha, names=FALSE)

    # Update the rejection matrix based on the 1-alpha critical value
    # Continue if there are any new rejections
    numrej = rejupdate(reject, diffmat, crit)
    newrej = numrej > 0

    # Display status
    cat("\r Step", j,". ", numrej, " new rejections.")
    utils::flush.console()
    j = j + 1
    gc()
  }
  return(reject)
}
