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

df = simfun(50, 500, 0.8, nullprop=0)

profvis::profvis({
  rankconf(df$yhat, df$sig2, type="BKFWER", k=500)
})

profvis::profvis({
  rankconf(df$yhat, df$sig2, type="BKFWER", k=500, thr=parallel::detectCores()-1)
})

test = list()
test[["Mogstad"]] = csranks::csranks(-df$yhat, sqrt(df$sig2))
test[["10-FWER"]] = rankconf(df$yhat, df$sig2, type="BKFWER", k=10, thr=parallel::detectCores()-1)
test[["500-FWER"]] = rankconf(df$yhat, df$sig2, type="BKFWER", k=100, thr=parallel::detectCores()-1)
test[["FDR"]] = rankconf(df$yhat, df$sig2, type="FDR", alpha=0.3)

ggplot(df, aes(x=rank(yhat), y=rank(y))) +
  geom_point() +
  geom_smooth(mapping=aes(x=rank(df$yhat), y=test[["Mogstad"]]$L), color="black") +
  geom_smooth(mapping=aes(x=rank(df$yhat), y=test[["Mogstad"]]$U), color="black") +
  geom_smooth(mapping=aes(x=rank(df$yhat), y=test[["500-FWER"]]$L), color="blue") +
  geom_smooth(mapping=aes(x=rank(df$yhat), y=test[["500-FWER"]]$U), color="blue") +
  geom_smooth(mapping=aes(x=rank(df$yhat), y=test[["FDR"]]$L, color="red")) +
  geom_smooth(mapping=aes(x=rank(df$yhat), y=test[["FDR"]]$U, color="red")) +
  theme_classic() +
  theme(
    axis.title.x = element_text("Estimated Rank"),
    axis.title.y = element_text("True Rank")
  )

sum(rank(df$y) >= test3$L & df$y <= test3$U)
# Sampfun ======================================================================
diffmat = matrix(selfouter(df$y, '-')/sqrt(selfouter(df$sig2, "+")), nrow(df), nrow(df))
diag(diffmat) = NA
sigmat = sqrt(selfouter(df$sig2, FUN = "+"))
ind = !is.na(diffmat)
sigmatind = sigmat[ind]

sampfun = function(sigmatind, ind, k, sd){
  booty = Rfast::Rnorm(length(sd), m=0, s=1)*sd
  return(
    kmax(
      abs(selfouter(booty, "-")[ind]/sigmatind),
      k
    )
  )
}

profvis::profvis({
  sampfun(sigmatind, ind, 5, sqrt(df$sig2))
})

# kmax =========================================================================
# Function to find the kth largest value of x
kmax1 = function(x, k){
  k = min(length(x), k)
  x[kit::topn(x, k, decreasing=T)[k]]
}
kmin1 = function(x, k){
  k = min(length(x), k)
  x[kit::topn(x, k, decreasing=F)[k]]
}
# Function to find the kth largest value of x
kmax2 = function(x, k){
  k = min(length(x), k)
  p = length(x) - k + 1
  return(
    sort(x, partial=p, decreasing=F)[p]
  )
}
kmin2 = function(x, k){
  k = min(length(x), k)
  p = length(x) - k + 1
  return(
    -sort(-x, partial=p, decreasing=F)[p]
  )
}

microbenchmark::microbenchmark(
  old = kmin1(x,300),
  new = kmin2(x,300)
)

# ==============================================================================
n = 20000
object.size(matrix(runif(n^2), n, n))


