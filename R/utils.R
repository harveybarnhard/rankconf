# Faster than using base::outer()
#' @export
selfouter = function(x, FUN="+"){
  FUN = match.fun(FUN)
  y = rep(x, rep.int(length(x), length(x)))
  FUN(x,y)
}
