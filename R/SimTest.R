as.SimTest <- 
function (RealValue, SimulatedValues) {
  st <- list(RealValue=RealValue, SimulatedValues=SimulatedValues)
  class(st) <- "SimTest"
  return(st)
}