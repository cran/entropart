plot.DivProfile <- 
function (x, ...) {

  # Save graphical parameters
  op <- par(no.readonly = TRUE)
	par (mfrow=c(2, 2))
	plot(y=x$TotalAlphaDiversity, x=x$Order, type="l", xlab="Order of Diversity", ylab="Alpha Diversity", main="Total Alpha Diversity", ...)
	plot(y=x$TotalBetaDiversity, x=x$Order, type="l", xlab="Order of Diversity", ylab="Beta Diversity", main="Beta Diversity", ...)
	plot(y=x$GammaDiversity, x=x$Order, type="l", xlab="Order of Diversity", ylab="Gamma Diversity", main="Gamma Diversity", ...)
	# Restore parameters
  par(op)
}
