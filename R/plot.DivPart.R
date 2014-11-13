plot.DivPart <- 
function (x, ...) {

  plot(c(0, x$GammaDiversity), c(0, length(x$CommunityAlphaDiversities)), type = "n", xlab = expression(paste(alpha, " and ", gamma, " diversity")), ylab = expression(paste(beta, " diversity")), ...)
  rect(0, 0, x$GammaDiversity, 1, lty=2)
  rect(0, 0, x$TotalAlphaDiversity, x$TotalBetaDiversity, lty=2)

}
