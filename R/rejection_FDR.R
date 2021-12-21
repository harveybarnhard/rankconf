# Benjamini-Yekutieli (FDR)
rej_FDR_BY = function(pvals, alpha, ...){
  n = nrow(pvals)
  m = n^2 - n
  cm = log(m) - digamma(1) + 1/(2*m)
  pranks = matrix(data.table::frank(c(pvals), ties.method="min", na.last=T), n, n)
  max_rank = max(c(pranks[pvals<=pranks*alpha/(m*cm)], -Inf), na.rm=T)
  return(pranks < max_rank)
}
