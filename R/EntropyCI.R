EntropyCI <-
function(FUN, Simulations = 100, Ns, ..., CheckArguments = TRUE) 
{
  if (CheckArguments) {
    CheckentropartArguments()
  }
  
  RealEst <- FUN(Ns, ..., CheckArguments = FALSE)
  
  # SimulateEntropy resamples the community and calculates FUN
  SimulateEntropy <- function(Progression) {
    # Draw Ns from a multinomial distribution
    SimNs <- rmultinom(1, sum(Ns), Ns)
    # FUN(simulated data)
    NewSim <- FUN(SimNs, ..., CheckArguments = FALSE)
    # update progress bar
    setTxtProgressBar(ProgressBar, Progression)
    return(NewSim)
  }
  
  # Simulate entropy
  ProgressBar <- txtProgressBar(min=0, max=Simulations)
  # Simulated values
  RawSimulatedEntropy <- sapply(1:Simulations, SimulateEntropy)
  # Recenter entropy
  return(RawSimulatedEntropy + RealEst - mean(RawSimulatedEntropy))
  
}
