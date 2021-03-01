
# TODO: Make this more efficient by just focusing on upper
#       triangular
rankconf = function(y,
                    sig2,
                    type="FDR",
                    method="BH",
                    alpha=0.05,
                    k=NA,
                    best="min"){
  # Handle Input ===============================================================
  if(any(is.na(y) | is.na(sig2) | sig2 <= 0)){
    stop("y and sig2 must not contain NA values, and sig2 must be positive")
  }

  # Make sure error types exist
  # TODO: Separate this between error rates and methods to control error rates
  type = toupper(type)
  method = toupper(method)
  if(!type%in%c("NAIVE", "BONF", "FDR", "FWER", "KFWER", "KBONF","PFER")){
    stop("method must be NAIVE, BONF, FDR, FWER, KFWER, PFER")
  }

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

  # Calculate naive p-values of all differences
  zscoremat = outer(y, y, '-')/sqrt(outer(sig2, sig2, "+"))
  pvals = 2*pnorm(-abs(zscoremat[upper.tri(zscoremat)]))

  # Find which tests to reject using the given method
  reject = do.call(paste0("rej", type), list(zscoremat, pvals, alpha, k))
  rejectmat = matrix(NA, nrow=n, ncol=n)
  rejectmat[upper.tri(rejectmat)] = reject
  rejectmat[lower.tri(rejectmat)] = t(rejectmat)[lower.tri(rejectmat)]
  diag(rejectmat) = FALSE

  # Create a "reject and positive difference" matrix.
  # A TRUE in the [i,j]th index indicates that the ith observation is
  # significantly better than the jth observation.
  diffpos = (zscoremat > 0) & rejectmat

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
rejNAIVE = function(zscoremat, pvals, alpha, k){
  return(pvals < alpha)
}

# Bonferroni (FWER)
rejBONF = function(zscoremat, pvals, alpha, k){
  return(pvals < alpha/length(pvals))
}

# Holm-Bonferroni (FWER)
rejFWER = function(scoremat, pvals, alpha, k){
  n = length(pvals)
  # Rank pvals from lowest to highest
  pranks = rank(pvals, ties.method="min")

  # Find first p value not low enough for rejection, reject everything below
  min_rank = min(c(pranks[pvals>=alpha/(n+1-pranks)], Inf))
  return(pranks < min_rank)
}

# Benjamini-Hochberg (FDR)
rejFDR = function(zscoremat, pvals, alpha, k){
  pranks = rank(pvals, ties.method="min")
  max_rank = max(c(pranks[pvals<=pranks*alpha/length(pvals)], -Inf))
  return(pranks < max_rank)
}

# Bonferroni (KBONF)
rejKBONF = function(zscoremat, pvals, alpha, k){
  return(pvals < k*alpha/length(pvals))
}

# Holm-Bonferroni (KFWER)
rejKFWER = function(zscoremat, pvals, alpha, k){
  n = length(pvals)
  # Rank pvals from lowest to highest
  pranks = rank(pvals, ties.method="min")

  # Create values against which to compare p values
  compvals = rep(NA, n)
  ind = pranks<=k
  compvals[ind]  = k*alpha/n
  compvals[!ind] = k*alpha/(n+k-pranks[!ind])

  # Find first p value not low enough for rejection, reject everything below
  min_rank = min(pranks[pvals>=compvals])
  return(pranks < min_rank)
}

# (REKFWER) https://arxiv.org/pdf/0710.2258.pdf Algorithm 2.1

# Function to find the kth largest value of x
kmax = function(x, k){
  n = length(x)
  sort(x, partial=n- k+ 1)[n-k+1]
}
sampfun = function(x, k){
  return(kmax(sample(x, size=length(x), replace=TRUE)))
}
# Function to find the bootstrap quantile of the kmax function over
# specific indices of x
pkmax = function(x, ind, quant, k, R=1000){
  # Set of indices
  kmaxdist = replicate(R, sampfun(x[rejind],k))
  return(quantile(kmaxdist, prob=quant))
}
rejBKFWER = function(zscoremat, pvals, alpha, k){
  tvals = abs(zscoremat[upper.tri(zscoremat)])
  reject = rep(0, length(tvals))

  # Step 1
  crit = pkmax(tvals, reject, 1-alpha, k)
  if(max(tvals)<=crit){
    return(reject > 0)
  }else{
    ind = tvals > crit
    newrej = TRUE
    reject[ind] = 1
  }
  if(sum(reject) < k){
    return(reject > 0)
  }

  # Find the k-1 least significant rejections
  kmax(tvals[reject], k-1)

  # Step 2,3,...
  j = 2
  while(newrej){
    crit = pkmax(tvals, recentrej | !reject, 1-alpha, k)
    ind = which(reject==0)
    reject[ind] = j*(tvals[ind] > crit)
    newrej = any(reject[ind]>0)
    j = j + 1
  }
  return(reject > 0)
}
