simfun = function(n, g, r, nullprop=0){
  sigvar = r/((1-r)*n)
  if(nullprop>0 & nullprop < 1){
    numnull = floor(nullprop*g)
    y = rnorm(g, mean=0, sd=sqrt(sigvar)/sqrt(1-nullprop))
    y[sample.int(g, numnull)] = 0
  }else if(nullprop==1){
    y = rep(0, g)
  }else{
    y = rnorm(g, mean=0, sd=sqrt(sigvar))
  }
  df = data.table::data.table(
    id = rep(1:g, 1, each=n),
    yi  = rnorm(n*g, mean=rep(y, 1, each=n), sd=1)
  )
  df = df[,
          list(yhat=mean(yi),
               sig2=var(yi)/n),
          by=id
  ]
  df[, id:=NULL]
  df$y = y
  return(df)
}

df = simfun(50, 400, 0.8, nullprop=0)

profvis::profvis({
  rankconf(df$yhat, df$sig2, type="BKFWER", k=500)
})

profvis::profvis({
  rankconf(df$yhat, df$sig2, type="BKFWER", k=500, thr=parallel::detectCores()-1)
})

test = list()
test[["Mogstad"]] = csranks::csranks(-df$yhat, sqrt(df$sig2))
test[["1-FWER"]] = rankconf(df$yhat, df$sig2, type="BKFWER", k=1, thr=parallel::detectCores()-1)
test[["10-FWER"]] = rankconf(df$yhat, df$sig2, type="BKFWER", k=10, thr=parallel::detectCores()-1)
test[["500-FWER"]] = rankconf(df$yhat, df$sig2, type="BKFWER", k=100, thr=parallel::detectCores()-1)
test[["FDR"]] = rankconf(df$yhat, df$sig2, type="FDR", alpha=0.05)

ggplot(df, aes(x=rank(yhat), y=rank(y))) +
  geom_point() +
  geom_point(mapping=aes(x=rank(df$yhat), y=test[["Mogstad"]]$L, color="black")) +
  geom_point(mapping=aes(x=rank(df$yhat), y=test[["Mogstad"]]$U, color="black")) +
  geom_point(mapping=aes(x=rank(df$yhat), y=test[["500-FWER"]]$L, color="blue")) +
  geom_point(mapping=aes(x=rank(df$yhat), y=test[["500-FWER"]]$U, color="blue")) +
  geom_point(mapping=aes(x=rank(df$yhat), y=test[["FDR"]]$L, color="red")) +
  geom_point(mapping=aes(x=rank(df$yhat), y=test[["FDR"]]$U, color="red")) +
  theme_classic() +
  xlab("Estimated Rank") + ylab("True Rank")

sum(rank(df$y) >= test3$L & df$y <= test3$U)

# Sampfun ======================================================================
diffmat = matrix(selfouter(df$y, '-')/sqrt(selfouter(df$sig2, "+")), nrow(df), nrow(df))
diag(diffmat) = NA
sigmat = sqrt(selfouter(df$sig2, FUN = "+"))
ind = sample(c(T, F), size=length(sigmat), replace=T, prob = c(0.8, 0.2))
sd = sqrt(df$sig2)
n  = nrow(df)

sampfun1 = function(sigmat, ind, k, sd){
  booty = Rfast::Rnorm(length(sd), m=0, s=1)*sd
  return(
    kmax(
      abs(selfouter(booty, "-")[ind]/sigmat[ind]),
      k
    )
  )
}
# Function that returns one bootstrap sample of the kth largest value
sampfun2 = function(sigmat, ind, k, distfun, ...){
  return(
    kmax(
      abs(selfouter(Rfast::Rnorm(length(sigmat), m=0, s=1)*sd, "-")[ind]/sigmat[ind]),
      k
    )
  )
}
# Function that returns one bootstrap sample of the kth largest value
sampfun3 = function(sigmat, ind, k, distfun, ...){
  return(
    kmax(
      abs(matrix(selfouter(Rfast::Rnorm(length(sd), m=0, s=1)*sd, "-"))/sigmat[ind]),
      k
    )
  )
}

microbenchmark::microbenchmark(
  sampfun1(sigmat, ind, 400, sd),
  sampfun2(sigmat, ind, 400, sd)
)

profvis::profvis({
  sampfun(sigmatind, ind, 800, sd)
})

# kmax =========================================================================
# k-largest element
kmax1 = function(x, k){
  k = min(length(x), k)
  if(k > 800){
    p = length(x) - k + 1
    sort(x, partial=p, decreasing=F)[p]
  }else{
    x[kit::topn(x, k, decreasing=T)[k]]
  }
}

# k-smallest element
kmin1 = function(x, k){
  k = min(length(x), k)
  if(k > 800){
    p = length(x) - k + 1
    -sort(-x, partial=p, decreasing=F)[p]
  }else{
    x[kit::topn(x, k, decreasing=F)[k]]
  }
}

# k-largest element
kmax3 = function(x, k){
  k = min(length(x), k)
  x[top_index(x, k)]
}

# k-smallest element
kmin3 = function(x, k){
  k = min(length(x), k)
  -x[top_index(-x, k)]
}

k=1000
microbenchmark::microbenchmark(
  oldmax = kmax1(df$y, k),
  newermax = kmax3(df$y, k),
  oldmin = kmin1(df$y, k),
  times = 1000
)


# ==============================================================================
nvec = c(75, 750, 7500)
res = rep(NA, length(nvec))
for(i in 1:length(nvec)){
  n = nvec[i]
  res[i] = object.size(matrix(runif(n^2), n, n))
}

# Sparse matrices ==============================================================
library(Matrix)
n=1000
m = list()
m[[1]] = matrix(sample(c(T,F), n^2, replace=T, prob=c(0.5,0.5)), n, n)
m[[2]] = Matrix(sample(c(T,F), n^2, replace=T, prob=c(0.5,0.5)), n, n)


# In-place matrix operations ===================================================

oldrej = function(c){
  logmat = matrix(c(F, F, F, F, F, F, F, F, F), 3, 3)
  diffmat = matrix(c(2, 1.2, 0.3, 1.2, 2, 0.7, 0.3, 0.7, 2), 3, 3)
  logmat[!logmat] = logmat[!logmat] > c
  diag(logmat) = F
}

newrej = function(c){
  logmat = matrix(c(F, F, F, F, F, F, F, F, F), 3, 3)
  diffmat = matrix(c(2, 1.2, 0.3, 1.2, 2, 0.7, 0.3, 0.7, 2), 3, 3)
  rejupdate(logmat, diffmat, c)
}

microbenchmark::microbenchmark(
  new = newrej(0.5),
  old = oldrej(0.5)
)
profmem(
  oldrej(0.5)
)
profmem(
  newrej(0.5)
)

# indupdate
oldind = function(n){
  logmat = matrix(rep(F, n^2), n, n)
  diffmat = matrix(runif(n^2), n, n)
  diffmat[lower.tri(diffmat)] = t(diffmat)[lower.tri(diffmat)]
  indmat = matrix(rep(F, n^2), n, n)
  indmat = (logmat & (indmat <= 0.7)) | !logmat
}
newind = function(n){
  logmat = matrix(rep(F, n^2), n, n)
  diffmat = matrix(runif(n^2), n, n)
  diffmat[lower.tri(diffmat)] = t(diffmat)[lower.tri(diffmat)]
  indmat = matrix(rep(F, n^2), n, n)
  indupdate(logmat, indmat, diffmat, 0.7, 3)
}

oldind2 = function(logmat, diffmat, indmat){
  indmat = (logmat & (indmat <= 0.7)) | !logmat
  return(indmat)
}
newind2 = function(logmat, diffmat, indmat){
  indupdate(logmat, indmat, diffmat, 0.7, 3)
}

n = 7500
logmat = matrix(rep(F, n^2), n, n)
diffmat = matrix(runif(n^2), n, n)
diffmat[lower.tri(diffmat)] = t(diffmat)[lower.tri(diffmat)]
indmat = matrix(rep(F, n^2), n, n)

# A little slower
microbenchmark::microbenchmark(
  new = newind(500),
  old = oldind(500)
)

# ...but negligible memory footprint
profmem(
  oldind2(logmat, diffmat, indmat)
)
profmem(
  newind(logmat, diffmat, indmat)
)






