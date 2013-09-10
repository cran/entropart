bcShannon<-
function(Ns, Correction="Best", CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  # No correction
  if (Correction == "None") {
    return(Shannon(Ns/sum(Ns), CheckArguments=FALSE))
  }
  
  # Eliminate 0
  Ns <- Ns[Ns > 0]
  Nrecords <- sum(Ns)
  SampleCoverage <- Coverage(Ns, CheckArguments=FALSE)
  if (Correction == "ChaoShen" | Correction == "Best") {
    CPs <- Ns/Nrecords*SampleCoverage 
    ChaoShen <- sum(-CPs*log(CPs)/(1-(1-CPs)^Nrecords))
  }
  if (Correction == "ChaoShen") {
    return(ChaoShen)  
  } 
  if (Correction == "Grassberger" | Correction == "Best") {
    # (-1)^n is problematic for long vectors (returns NA for large values). It is replaced by 1-n%%2*2 (Ns is rounded if is not an integer)
    Grassberger <- sum(Ns/Nrecords*(log(Nrecords)-digamma(Ns)-(1-round(Ns)%%2*2)/(Ns+1)))
  }
  if (Correction == "Grassberger") {
    return(Grassberger)
  }
  if (Correction == "Best") {
    return(max(ChaoShen, Grassberger))
  }
  if (Correction == "Grassberger2003" | Correction == "Schurmann") {
    # Define a function to calculate the integral in the bias formuma for each value of N
    Integral <- function(n, upper) integrate(function(t, n) t^(n-1)/(1+t), 0, upper, n) 
  }
  if (Correction == "Grassberger2003") {
    Integral.V <- unlist(sapply(Ns, Integral, upper = 1)["value",])
  }
  if (Correction == "Schurmann") {
    Integral.V <- unlist(sapply(Ns, Integral, upper = exp(-1/2))["value",])
  }
  if (Correction == "Grassberger2003" | Correction == "Schurmann") {
    return(sum(Ns/Nrecords*(digamma(Nrecords)-digamma(Ns)-(1-Ns%%2*2)*Integral.V)))
  }
  if (Correction == "Holste" | Correction == "Bonachela") {
    seql <- 1:(length(Ns)+Nrecords)
    invl <- 1/seql
    cumul <- function(n) {sum(invl[n:length(invl)])}
    suminvl <- sapply(seql, cumul) 
  if (Correction == "Holste") {
    return(sum((Ns+1)/(length(Ns)+Nrecords)*suminvl[Ns+2])) 
   } else {
    return(sum((Ns+1)/(2+Nrecords)*suminvl[Ns+2]))    
   }
  } 
  warning("Correction was not recognized")
  return (NA)
}
