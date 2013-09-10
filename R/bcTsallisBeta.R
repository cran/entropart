bcTsallisBeta <-
function(Ns, Nexp = NULL, q, Correction = "Best", CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  if (length(Ns) != length(Nexp)) {
    stop("Ns and Nexp should have the same length.")
  }  
    
  # No correction
  if (Correction == "None") {
    return(TsallisBeta(Ns/sum(Ns), Nexp/sum(Nexp), q, CheckArguments=FALSE))
  }
  
  # Sample coverage
  Nrecords <- sum(Ns)
  SampleCoverage <- Coverage(Ns, CheckArguments=FALSE)
  # Sample coverage (expected)
  Nrecordsexp <- sum(Nexp)
  SampleCoverageexp <- Coverage(Nexp, CheckArguments=FALSE)
  if (Correction == "ChaoShen" | Correction == "Best") {
    CiPsi <- SampleCoverage * Ns / Nrecords
    CPs <- SampleCoverageexp * Nexp / Nrecordsexp
    dataBeta <- CiPsi^q * lnq(CiPsi/CPs, q) / (1 -(1-CiPsi)^Nrecords)    
    # force 0log0=0                                                                                                                                       
    dataBeta[Ns == 0] <- 0
    return (sum(dataBeta))  
  } 
  warning("Correction was not recognized")
  return (NA)
}
