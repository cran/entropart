bcHqzBeta <-
function(Ns, Nexp = NULL, q = 1, Z = diag(length(Ns)), Correction = "Best", CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  if (length(Ns) != length(Nexp)) {
    stop("Ns and Nexp should have the same length.")
  }  
    
  # No correction available yet
  if (Correction == "None" | Correction == "Best") {
    return(HqzBeta(Ns/sum(Ns), Nexp/sum(Nexp), q, Z, CheckArguments=FALSE))
  }

  warning("Correction was not recognized")
  return (NA)
}
