DivEst <-
function(q = 0, MC, Biased = TRUE, Correction = "Best", Tree = NULL, Normalize = TRUE, Z = NULL, Simulations = 100, CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  # Preprocess the tree
  ppTree <- Preprocess.Tree(Tree)
  if (Normalize) {
    Height <- 1
  } else {
    Height <- ppTree$Height
  }  

  # Estimation from data
  RealEst <- DivPart(q, MC, Biased, Correction, ppTree, Normalize, Z, CheckArguments=FALSE)

  # RedrawSpecies resamples a community according to species abundances.
  RedrawSpecies<- function(SpeciesAbundances){
    rmultinom(1, sum(SpeciesAbundances), SpeciesAbundances)
  }

  # SimulateEntropy resamples all communities and calculates entropy
  SimulateEntropy <- function(Progression) {
    SimNsi <- apply(MC$Nsi, 2, RedrawSpecies)
    rownames(SimNsi) <- rownames(MC$Nsi) 
    SimMC <- Preprocess.MC(SimNsi, MC$Wi)
    NewSim <- DivPart(q, SimMC, Biased, Correction, Tree, Normalize, Z, CheckArguments=FALSE)
    # update progress bar
    setTxtProgressBar(ProgressBar, Progression)
    c(NewSim$TotalAlphaEntropy, NewSim$TotalBetaEntropy, NewSim$GammaEntropy)
  }
  
  # Simulate entropy
  ProgressBar <- txtProgressBar(min=0, max=Simulations)
  RawSimulatedEntropy <- sapply(1:Simulations, SimulateEntropy)
  # Recenter entropy
  SimulatedEntropy <- RawSimulatedEntropy+c(RealEst$TotalAlphaEntropy, RealEst$TotalBetaEntropy, RealEst$GammaEntropy)-apply(RawSimulatedEntropy ,1, mean)
  # Transform entropy to diversity
  SimulatedDiversity <- SimulatedEntropy
  if (q == 1) { 
    SimulatedDiversity <- exp(SimulatedEntropy)
  } else {
    SimulatedDiversity[1,] <- expq(SimulatedEntropy[1,] / Height, q) * Height
    SimulatedDiversity[2,] <- expq(SimulatedEntropy[2,] / Height / (1 - (q-1)*SimulatedEntropy[2,]/Height), q) * Height
    SimulatedDiversity[3,] <- expq(SimulatedEntropy[3,] / Height, q) * Height
  }
  dimnames(SimulatedDiversity)[[1]] <- list("Alpha", "Beta", "Gamma")

  DivEst <- RealEst
  DivEst$SimulatedDiversity <- SimulatedDiversity
  class(DivEst) <- c("DivEst", class(RealEst))

  return (DivEst)
}
