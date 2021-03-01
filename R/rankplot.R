rankplot = function(y, L, U) {
  ranky = rank(y, ties.method="min")
  g = ggplot(data=data.frame(ranky, L, U)) +
    geom_linerange(mapping=aes(x=ranky, ymin=L, ymax=U)) +
    theme_classic()
  return(g)
}
