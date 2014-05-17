plot.DivProfile <- 
function (x, ..., main = NULL, xlab = "Order of Diversity", ylab = NULL, Which = "All") {

  # Save graphical parameters
  if (Which == "All") {
    op <- par(no.readonly = TRUE)
    par (mfrow=c(2, 2))    
  }
  if (Which == "All" | (Which == "Alpha" & is.null(main))) main <- "Total Alpha Diversity"
  if (Which == "All" | (Which == "Alpha" & is.null(ylab))) ylab <- expression(paste(alpha, " diversity"))
  if (Which == "All" | Which == "Alpha") {
    plot(y=x$TotalAlphaDiversity, x=x$Order, type="l", main=main, xlab=xlab, ylab=ylab, ...)
  }
  if (Which == "All" | (Which == "Communities" & is.null(main))) main <- "Alpha Diversity of Communities"
  if (Which == "All" | (Which == "Communities" & is.null(ylab))) ylab <- expression(paste(alpha, " diversity"))
  if (Which == "All" | Which == "Communities") {
    plot(x$CommunityAlphaDiversities[, 1] ~ x$Order, type="n", xlim=c(min(x$Order), max(x$Order)), ylim=c(min(x$CommunityAlphaDiversities), max(x$CommunityAlphaDiversities)), main=main, xlab=xlab, ylab=ylab, ...)
    for (Community in (1:ncol(x$CommunityAlphaDiversities))) {
      lines(x=x$Order, y=x$CommunityAlphaDiversities[, Community], lty=Community)
    }  
    if (Which == "Communities") {
      legend("topright", colnames(x$CommunityAlphaDiversities), lty=1:ncol(x$CommunityAlphaDiversities), inset=0.01)
    }
  }
  if (Which == "All" | (Which == "Beta" & is.null(main))) main <- "Beta Diversity"
  if (Which == "All" | (Which == "Beta" & is.null(ylab))) ylab <- expression(paste(beta, " diversity"))
  if (Which == "All" | Which == "Beta") {
    plot(y=x$TotalBetaDiversity, x=x$Order, type="l", main=main, xlab=xlab, ylab=ylab, ...)
  }
  if (Which == "All" | (Which == "Gamma" & is.null(main))) main <- "Gamma Diversity"
  if (Which == "All" | (Which == "Gamma" & is.null(ylab))) ylab <- expression(paste(gamma, " diversity"))
  if (Which == "All" | Which == "Gamma") {
    plot(y=x$GammaDiversity, x=x$Order, type="l", main=main, xlab=xlab, ylab=ylab, ...)
  }
	# Restore parameters
  if (Which == "All") {
    par(op)
  }
}
