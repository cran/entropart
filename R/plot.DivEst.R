plot.DivEst <- 
function (x, ..., main = NULL, Which = "All") {

	# Save graphical parameters
	op <- par(no.readonly = TRUE)
  par (mfrow=c(2, 2))
  
	if (Which == "All" | (Which == "Alpha" & is.null(main))) main <- "Alpha Diversity"
	if (Which == "All" | Which == "Alpha") {
	  plot(as.SimTest(x$TotalAlphaDiversity, x$SimulatedDiversity["Alpha",]), main=main, ...)
	}
	if (Which == "All" | (Which == "Beta" & is.null(main))) main <- "Beta Diversity"
	if (Which == "All" | Which == "Beta") {
	  plot(as.SimTest(x$TotalBetaDiversity, x$SimulatedDiversity["Beta",]), main=main, ...)
	}
	if (Which == "All" | (Which == "Gamma" & is.null(main))) main <- "Gamma Diversity"
	if (Which == "All" | Which == "Gamma") {
	  plot(as.SimTest(x$GammaDiversity, x$SimulatedDiversity["Gamma",]), main=main, ...)
	}
	
	# Legend and restore parameters
	if (Which == "All") {
    par(mar=c(0, 0, 0, 0))
	  plot(0:10, 0:10, type="n", xlab=NULL, frame.plot=FALSE, xaxt="n", yaxt="n", col.lab="white")
	  leg <- c("Null Distribution", "True Estimate", "95% confidence interval") 
	  legend(2, 8, leg, col = c(1, 2, 1), lty = 1:3, merge = TRUE, cex=1)
	  par(op)
	}
}
