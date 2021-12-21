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
#' @param thr Number of threads to use for parallel processing. Number of cores minus one by default.
#' @export
rankconf = function(y,
                    sig2,
                    type="PCER",
                    method="NAIVE",
                    alpha=0.05,
                    k=1,
                    best="min",
                    thr=parallel::detectCores()-1
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

  # Calculate naive one-sided p-values of all differences
  diffmat = matrix(selfouter(y, '-')/sqrt(selfouter(sig2, "+")), n, n)
  diag(diffmat) = NA
  if(!(type=="FWER" & method=="R")){
    pvals = 1-pnorm(diffmat)
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
  lowerrank = rowSums(diffmat) + 1
  upperrank = n - colSums(diffmat)
  return(list(
    L = lowerrank,
    U = upperrank
  ))
}
