Coverage <-
function(Ns, Estimator = "ZhangHuang", CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  # Round values
  Ns <- as.integer(Ns)
  # Eliminate zeros
  Ns <- Ns[Ns>0]
  # Calculate abundance distribution
  DistN <- tapply(Ns, Ns, length)
  Singletons <- DistN["1"]

  # No singletons, C=1
  if (is.na(Singletons)) {
    return(1)
  }
  
  # Abundances
  Nu <- as.integer(names(DistN))
  SampleSize <- sum(Ns)

  # Singletons only
  if (Singletons == SampleSize) {
    warning ("Sample coverage is 0, most bias corrections will return NaN.")
    return(0)
  }

  if (Estimator == "ZhangHuang") {
    Ps <- Ns/SampleSize
    if (any(Ps >= .5)) {
      warning ("Zhang-Huang sample coverage cannot be estimated because one probability is over 1/2. Chao estimator is returned.")
      Estimator <- "Chao"
    } else {
      # Use Nu%%2*2-1 for (-1)^(Nu+1)
      return(1 - sum((Nu%%2*2-1) / choose(SampleSize, Nu) * DistN))
    }    
  }
  if (Estimator == "Chao") {
    Doubletons <- DistN["2"]
    if (is.na(Doubletons)) {
      warning ("Chao's sample coverage cannot be estimated because there are no doubletons. Turing estimator is returned.")
      Estimator <- "Turing"
    } else {
    return(1 - Singletons / SampleSize *((SampleSize - 1) * Singletons / ((SampleSize - 1) * Singletons + 2 * Doubletons)))
    }
  }
  if (Estimator == "Turing") {
    return(1 - Singletons / SampleSize)
  }
  warning("Correction was not recognized")
  return (NA)

}
