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
                    thr=parallel::detectCores()-1){
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

kmax = function(x, k){
  k = min(length(x), k)
  if(k > 800){
    p = length(x) - k + 1
    sort(x, partial=p, decreasing=F)[p]
  }else{
    x[kit::topn(x, k, decreasing=T)[k]]
  }
}
kmin = function(x, k){
  k = min(length(x), k)
  if(k > 800){
    p = length(x) - k + 1
    -sort(-x, partial=p, decreasing=F)[p]
  }else{
    x[kit::topn(x, k, decreasing=F)[k]]
  }
}

# Function that returns one bootstrap sample of the kth largest value
sampfun = function(sigmatind, ind, k, distfun, ...){
  booty = do.call(distfun, list(...))
  return(
    kmax(
      abs(selfouter(booty, "-")[ind]/sigmatind),
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
  on.exit(parallel::stopCluster(cl))

  # Initialize with no rejections
  reject = matrix(0, nrow=length(sig2), ncol=length(sig2))

  # Test statistics and SDs
  diffmat = abs(diffmat)
  diag(diffmat) = NA
  sigmat = sqrt(selfouter(sig2, FUN = "+"))

  # Step 1, calculate 1-alpha critical values
  #     a) Reject all values greater than critical value
  #     b) If k or fewer observations are rejected, then stop
  ind = !is.na(diffmat)
  n = nrow(sigmat)
  s = sqrt(sig2)
  parallel::clusterExport(
    cl, varlist=c("sigmat","ind", "k", "sampfun", "kmax",
                  "n", "s","distfun", "selfouter"),
    envir=environment()
  )
  kmaxdist = parallel::parSapply(
    cl=cl,
    X=rep(list(1), R),
    FUN=function(x) do.call(
      sampfun, list(sigmatind=sigmat[ind], ind=ind, k=k, distfun=distfun, sd=s)
    )
  )
  # Return the quantile of the distribution
  crit = quantile(kmaxdist, probs=1-alpha, names=FALSE)
  newrej = TRUE
  reject = diffmat > crit
  diag(reject) = FALSE

  if(sum(reject, na.rm=TRUE) < k){
    return(reject > 0)
  }
  cat("\r Step 1 Total Rejected", sum(reject>0, na.rm=T))

  # Step 2,3,...
  j = 2
  while(newrej & k>1){
    # Most "recent" k-1 rejections. Equivalently, the k-1 least significant
    # tests that have already been rejected
    if(k>1){
      recentrej = !reject
    }else{
      recentrej = reject & (diffmat <= kmin(x=diffmat[reject], k=k-1))
    }

    # Calculate critical value over unrejected tests and k-1 least significant
    # tests that have already been rejected
    ind = (!reject | recentrej) & !is.na(diffmat)
    parallel::clusterExport(
      cl, varlist="ind",
      envir=environment()
    )
    kmaxdist = parallel::parSapply(
      cl=cl,
      X=rep(list(1), R),
      FUN=function(x) do.call(
        sampfun, list(sigmatind=sigmat[ind], ind=ind, k=k, distfun=distfun, sd=s)
      )
    )
    # Return the quantile of the distribution
    crit = quantile(kmaxdist, probs=1-alpha, names=FALSE)
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
