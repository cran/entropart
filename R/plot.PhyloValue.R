plot.PhyloValue <- 
function (x, xlab = "T", ylab = NULL, main = NULL, ...) {
  
  Entity <- ""
  # Entity
  if (is.PhyloEntropy(x)) {
    Entity <- "Entropy"
  } else {
    if (is.PhyloDiversity(x)) {
      Entity <- "Diversity"
    }
  }

  # ylab
  if (is.null(ylab))
    ylab <- Entity
  
  # main
  if (is.null(main))
    main <- paste(Entity, "along the tree")
  
  plot(x=c(0, names(x$Cuts)), y=c(x$Cuts, 0), type="s", xlab=xlab, ylab=ylab, main=main, ...)
  abline(h=x$Total, lty=2)
  
}
