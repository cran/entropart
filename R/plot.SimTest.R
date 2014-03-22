plot.SimTest <- 
function (x, Quantiles = c(0.025, 0.975), ..., colValue = "red", lwdValue = 2, ltyValue = 2, colQuantiles = "black", lwdQuantiles = 1, ltyQuantiles = 2)
{
  plot(density(x$SimulatedValues), ...)
  abline(v=x$RealValue, col=colValue, lwd=lwdValue, lty=ltyValue)
  for (qt in Quantiles) {
    abline(v=quantile(x$SimulatedValues, probs = qt), col=colQuantiles, lwd=lwdQuantiles, lty=ltyQuantiles)
  }
}
