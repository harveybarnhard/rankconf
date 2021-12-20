# Naive method (PCER)
rej_PCER_NAIVE = function(pvals, alpha, ...){
  return(pvals < alpha/2)
}
