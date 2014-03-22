summary.SimTest <-
function(object, Quantiles = c(0.025, 0.975), ...) {
  
  cat("Real value: ", object$RealValue, "\n")
  cat("Quantile in the simulated distribution: ", ecdf(object$SimulatedValues)(object$RealValue), "\n")
  
  cat("Quantiles of simulations:\n")
  for (qt in Quantiles) {
    cat(sprintf("%1.2f%%", 100*qt), ": ", quantile(object$SimulatedValues, probs = qt), "\n")
  }
  cat("Mean simulated value: ", mean(object$SimulatedValues))
  return(invisible(NULL))
}