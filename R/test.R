x = rnorm(400, mean=0, sd=8)
sig2 = rnorm(400, mean=4, sd=2)^2
numfun = function(lower, upper){
  gapsize = upper-lower
  results = hb_rank(x, sig2, 0.8, 10*gapsize)
  numbelow = sum(results$U < lower)
  numabove = sum(results$L > upper)
  return(sum(c(numbelow, numabove)))
}

n=198
result = rep(NA, n)
for(i in 1:n){
  result[i] = numfun(i, length(x)-i)
  cat(i, "\n")
}
plot(result)
