plot.MCentropy <- 
function (x, ...) {

  barplot(c(x$Communities, 0, x$Total),
          beside = TRUE,
          width = c(x$Weights, .5, 1),
          names.arg = c(names(x$Communities), "", "Metacommunity"),
          ylab = "Entropy",
          ...
  )

}
