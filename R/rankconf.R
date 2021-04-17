
# TODO: Make this more efficient by just focusing on upper
#       triangular
rankconf = function(y,
                    sig2,
                    type="FDR",
                    alpha=0.05,
                    k=NA,
                    best="min",
                    thr=1){
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
    message("Rank 1 corresponds to the lowest estimate")
  }else if(best=="max"){
    message("Rank 1 corresponds to the highest estimate")
    y = -y
  }

  # Calculate naive one-sided p-values of all differences
  diffmat = outer(y, y, '-')/sqrt(outer(sig2, sig2, "+"))
  pvals = 1-pnorm(diffmat)

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
  diffpos = (diffmat > 0) & reject & !is.na(diffmat)

  # Calculate number of ranks above and below.
  # Summing over the jth column is equivalent to adding up the number of
  # observations that are significantly worse than the jth observation.
  # Summing over the ith row is equivalent to adding up the number of
  # observations that are significantly better than the ith observation.
  lowerrank = rowSums(diffpos) + 1
  upperrank = n - colSums(diffpos)
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
  pranks = matrix(rank(pvals, ties.method="min", na.last=T), n, n)

  # Find first p value not low enough for rejection, reject everything below
  min_rank = min(c(pranks[pvals>=alpha/(m +1-pranks)], Inf), na.rm=T)
  return(pranks < min_rank)
}

# Benjamini-Yekutieli (FDR)
rejFDR = function(pvals, alpha, ...){
  n = nrow(pvals)
  m = n^2 - n
  cm = log(m) -digamma(1) + 1/(2*m)
  pranks = matrix(rank(pvals, ties.method="min", na.last=T), n, n)
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
  pranks = matrix(rank(pvals, ties.method="min", na.last=T), n, n)

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

# Function to find the kth largest value of x
kmax = function(x, k){
  k = min(length(x), k)
  return(
    x[kit::topn(x, k)[k]]
  )
}
kmin = function(x, k){
  k = min(length(x), k)
  return(
    x[kit::topn(x, k, decreasing=FALSE)[k]]
  )
}

# Function that returns one bootstrap sample
sampfun = function(sigmat, ind, k, distfun, ...){
  booty = do.call(distfun, list(...))
  return(
    kmax(
      abs(outer(booty, booty, "-")[ind]/sigmat[ind]),
      k
    )
  )
}
# Function to find the bootstrap quantile of the kmax function over
# specific indices of x
pkmax = function(sigmat, ind, quant, k, R, distfun, thr, ...) with(list(...), {
  # Make R bootstrap samples, calculating the kth largest
  # observation in each sample
  cl = parallel::makeCluster(thr)
  parallel::clusterExport(
    cl, varlist=c("sigmat","ind", "quant", "k", "R", "sampfun", "kmax",
                  names(list(...))),
    envir=environment()
  )
  kmaxdist = parallel::parSapply(
    cl=cl,
    X=rep(list(1), R),
    FUN=function(x) do.call(
      sampfun,
      c(
        list(sigmat=sigmat, ind=ind, k=k, distfun=distfun),
        list(...)
      )
    )
  )
  parallel::stopCluster(cl)

  # Return the quantiile of the distrbution
  return(quantile(kmaxdist, probs=quant, names=FALSE))
})

# Main function
rejBKFWER = function(diffmat, sig2, alpha, k, R=1000, distfun="rnorm", thr=0, ...){
  # Initialize with no rejections
  reject = matrix(0, nrow=length(sig2), ncol=length(sig2))

  # Test statistics and SDs
  diffmat = abs(diffmat)
  diag(diffmat) = NA
  sigmat = sqrt(outer(sig2, sig2, FUN = "+"))

  # Step 1, calculate 1-alpha critical values
  #     a) Reject all values greater than critical value
  #     b) If k or fewer observations are rejected, then stop
  crit = pkmax(
    sigmat = sigmat,
    ind    = !is.na(diffmat),
    quant  = 1-alpha,
    k = k, R=R, thr=thr,
    distfun = distfun, n = nrow(sigmat), sd = sqrt(sig2)
  )
  newrej = TRUE
  reject = diffmat > crit
  diag(reject) = FALSE

  if(sum(reject, na.rm=TRUE) < k){
    return(reject > 0)
  }
  cat("\r Step 1 Total Rejected", sum(reject>0, na.rm=T))

  # Step 2,3,...
  j = 2
  while(newrej){
    # Most "recent" k-1 rejections. Equivalently, the k-1 least significant
    # tests that have already been rejected
    if(k>1){
      recentrej = !reject
    }else{
      recentrej = reject & (diffmat <= kmin(x=diffmat[reject], k=k-1))
    }

    # Calculate critical value over unrejected tests and k-1 least significant
    # tests that have already been rejected
    crit = pkmax(
      sigmat = sigmat,
      ind    = (!reject | recentrej) & !is.na(diffmat),
      quant  = 1-alpha,
      k = k, R=R, thr=thr,
      distfun = distfun, n = nrow(sigmat), sd = sqrt(sig2)
    )
    reject[!reject] = j*(diffmat[!reject] > crit)
    diag(reject) = FALSE
    if(!any(reject==j)){
      newrej = FALSE
    }
    cat("\r Step",j, "Total Rejected", sum(reject>0, na.rm=T))
    flush.console()
    j = j + 1
  }
  return(reject > 0)
}
