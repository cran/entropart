MetaCommunity <-
function(Abundances, Weights)
{
  Nspecies <- length(Abundances[, 1])
  if (is.factor(Abundances[,1])) {
    FirstColumnOfData <- 2
    SpeciesNames <- Abundances[,1]
    Ncommunities <- length(Abundances[1, ])-1
  } else {
    FirstColumnOfData <- 1
    SpeciesNames <- as.factor(rownames(Abundances))
    Ncommunities <- length(Abundances[1, ])
  }
  # Matrix containing p_si
  Nsi <- as.matrix(Abundances[, FirstColumnOfData:length(Abundances[1, ])])
  dimnames(Nsi)[[1]] <- SpeciesNames
  # Vector of weights
  if (is.vector(Weights)) {
    Wi <- Weights/sum(Weights)
  } else {
    Wi <- Weights$Weights/sum(Weights$Weights)
  }
  # Name the weight vector
  names(Wi) <- colnames(Nsi)
  Preprocess.MC(Nsi, Wi)
}
