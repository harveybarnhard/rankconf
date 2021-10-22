rankplot = function(y, L, U, size=1, title="", highlight=c(0,0)) {
  ranky = rank(y, ties.method="min")
  uptail  = length(y) - highlight[2]*length(y)
  lowtail = highlight[1]*length(y)
  color = rep("Neither", length(y))
  color[U < lowtail] = "good"
  color[L > uptail]  = "bad"
  g = ggplot2::ggplot(data=data.frame(ranky, L, U, color)) +
    ggplot2::geom_linerange(mapping=aes(x=ranky, ymin=L, ymax=U, colour=color),
                   size=size) +
    ggplot2::xlab("") + ggplot2::ylab("") +
    ggplot2::scale_x_continuous(expand=c(0.01,0.01)) +
    ggplot2::scale_y_continuous(expand=c(0.01,0.01)) +
    ggplot2::scale_color_manual(values=c("#1a4c8b", "#E12900", "grey"), guide = FALSE) +
    ggplot2::theme_classic() +
    ggplot2::ggtitle(label=title) +
    ggplot2::theme(
      text = ggplot2::element_text(family="LM Roman"),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
  return(g)
}

