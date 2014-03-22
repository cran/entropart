MetaCommunity <-
function(Abundances, Weights = rep(1, ncol(Abundances)))
{
  Nspecies <- length(Abundances[, 1])
  if (is.factor(Abundances[,1])) {
    FirstColumnOfData <- 2
    SpeciesNames <- Abundances[,1]
    Ncommunities <- length(Abundances[1, ])-1
  } else {
    FirstColumnOfData <- 1
    if (is.null(rownames(Abundances))) {
      # Create species names
      SpeciesNames <- as.factor(paste("sp", 1:(nrow(Abundances)), sep=""))
    } else {
      # Read species names
      SpeciesNames <- as.factor(rownames(Abundances))
    }
    Ncommunities <- length(Abundances[1, ])
  }
  
  # Community names
  if (is.null(colnames(Abundances))) {
    # Create community names
    colnames(Abundances) <- paste("P", 1:(ncol(Abundances)), sep="")
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
