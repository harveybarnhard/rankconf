#' kthe smallest element
#' @description Finds the kth smallest element of a numeric vector.
#' @keywords internal
kmin = function(x, k){
  k = min(length(x), k)
  if(k > 800){
    p = length(x) - k + 1
    -sort(-x, partial=p, decreasing=F)[p]
  }else{
    x[kit::topn(x, k, decreasing=F)[k]]
  }
}



#' ~3x as fast relative to base::outer()
#' @description Produces an outer product of a single vector with itself.
#' @keywords internal
selfouter = function(x, FUN="+"){
  match.fun(FUN)(x, rep(x, rep.int(length(x), length(x))))
}
