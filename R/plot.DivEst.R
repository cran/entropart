plot.DivEst <- 
function (x, ...) {

	# Save graphical parameters
	op <- par(no.readonly = TRUE)
  par (mfrow=c(2, 2))
  
  plot(density(x$SimulatedDiversity["Alpha",]), col="black", lwd=2, main="Alpha Diversity", xlab ="Diversity", ...)
  abline(v=x$TotalAlphaDiversity, col="red", lwd=2, lty=2)
  abline(v=quantile(x$SimulatedDiversity["Alpha",], probs = 0.025), col="black", lwd=1, lty=3)
  abline(v=quantile(x$SimulatedDiversity["Alpha",], probs = 0.975), col="black", lwd=1, lty=3)
  
  plot(density(x$SimulatedDiversity["Beta",]), col="black", lwd=2, main="Beta Diversity", xlab ="Diversity", ...)
  abline(v=x$TotalBetaDiversity, col="red", lwd=2, lty=2)
  abline(v=quantile(x$SimulatedDiversity["Beta",], probs = 0.025), col="black", lwd=1, lty=3)
  abline(v=quantile(x$SimulatedDiversity["Beta",], probs = 0.975), col="black", lwd=1, lty=3)

  plot(density(x$SimulatedDiversity["Gamma",]), col="black", lwd=2, main="Gamma Diversity", xlab ="Diversity", ...)
  abline(v=x$GammaDiversity, col="red", lwd=2, lty=2)
  abline(v=quantile(x$SimulatedDiversity["Gamma",], probs = 0.025), col="black", lwd=1, lty=3)
  abline(v=quantile(x$SimulatedDiversity["Gamma",], probs = 0.975), col="black", lwd=1, lty=3)
  
  plot(0:10, 0:10, type="n", xlab=NULL, frame.plot=FALSE, xaxt="n", yaxt="n", col.lab="white")
  leg <- c("Null Distribution", "True Estimate", "95% confidence interval") 
  legend(0, 10, leg, col = c(1, 2, 1), lty = 1:3, merge = TRUE, cex=1)
	
  # Restore parameters
	par(op)
}
