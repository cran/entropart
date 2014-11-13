bcShannon<-
function(Ns, Correction = "Best", CheckArguments = TRUE) 
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
  if (Correction == "ChaoShen") {
    SampleCoverage <- Coverage(Ns, CheckArguments=FALSE)
    CPs <- Ns/Nrecords*SampleCoverage 
    ChaoShen <- sum(-CPs*log(CPs)/(1-(1-CPs)^Nrecords))
  }
  if (Correction == "ChaoShen") {
    return(ChaoShen)  
  } 
  if (Correction == "Grassberger") {
    # (-1)^n is problematic for long vectors (returns NA for large values). It is replaced by 1-n%%2*2 (Ns is rounded if is not an integer)
    Grassberger <- sum(Ns/Nrecords*(log(Nrecords)-digamma(Ns)-(1-round(Ns)%%2*2)/(Ns+1)))
  }
  if (Correction == "Grassberger") {
    return(Grassberger)
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
  if (Correction == "ChaoWangJost" | Correction == "Best") {
    # Calculate abundance distribution
    DistN <- tapply(Ns, Ns, length)
    Singletons <- DistN["1"]
    Doubletons <- DistN["2"]
    # Calculate A
    if (is.na(Doubletons)) {
      if (is.na(Singletons)) {
        A <- 1
      } else {
        A <- 2/((Nrecords-1)*(Singletons-1)+2)
      }
    } else {
      A <- 2*Doubletons/((Nrecords-1)*Singletons+2*Doubletons)
    }
    Part2 <- sapply(1:(Nrecords-1), function(r) 1/r*(1-A)^r) 
    ChaoWangJost <- sum(Ns/Nrecords*(digamma(Nrecords)-digamma(Ns)))
    if (!is.na(Singletons) & (A != 1)) {
      ChaoWangJost <- as.numeric(ChaoWangJost + Singletons/Nrecords*(1-A)^(1-Nrecords)*(-log(A)-sum(Part2)))
    }
    return(ChaoWangJost)  
  }
  if (Correction == "ZhangHz") {
    # Values of v
    V <- 1:(Nrecords-1)
    Ps <- Ns/Nrecords
    # Weight part. Taken in log or goes to Inf for v > 1000; gamma cannot be used for large n, lgamma is preferred.
    lnw_v <- ((V+1)*log(Nrecords)+lgamma(Nrecords-V)-lgamma(Nrecords+1)-log(V))
    # p_V_Ps is an array, containing (1 - p_s - j/n) for each species (lines) and all j from 0 to n-2. Because array indexation starts from 1 in R, j is replaced by j-1.
    p_V_Ps <- outer(Ps, V, function(Ps, j) 1 -Ps -(j-1)/Nrecords)
    # Useful values are products from j=0 to v, so prepare cumulative products
    p_V_Ps <- t(apply(p_V_Ps, 1, cumprod))
    # Sum of products, weighted by p_s
    S_s <- function(v) {
      sum(Ps*p_V_Ps[1:length(Ps), v])
    }
    # Apply S_s to all values of v. Use logs or w_v goes to Inf.
    return(sum(exp(lnw_v + log(sapply(V, S_s)))))
  }
  warning("Correction was not recognized")
  return (NA)
}
