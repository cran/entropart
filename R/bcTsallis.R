bcTsallis <-
function(Ns, q = 1, Correction = "Best", CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  # No correction
  if (Correction == "None") {
    return(Tsallis(Ns/sum(Ns), q, CheckArguments = FALSE))
  }

  # Eliminate 0
  Ns <- Ns[Ns > 0]
  N <- sum(Ns)

  # Common code for ZhangGrabchak
  if (Correction == "ZhangGrabchak") {
    Ps <- Ns/N
    V <- 1:(N-1)
    # p_V_Ns is an array, containing (1 - (n_s-1)/(n-j)) for each species (lines) and all j from 1 to n-1
    p_V_Ns <- outer(Ns, V, function(Ns, j) 1- (Ns-1)/(N-j))
    # Useful values are products from j=1 to v, so prepare cumulative products
    p_V_Ns <- t(apply(p_V_Ns, 1, cumprod))
    # Sum of products weighted by w_v
    S_v <- function(s) {
      Usedv <- 1:(N-Ns[s])
      return(sum(w_v[Usedv]*p_V_Ns[s, Usedv]))
    }
  }
  
  # Shannon
  if (q == 1) {
    if (Correction == "Best") {
      ChaoShen <- bcShannon(Ns, Correction="ChaoShen", CheckArguments=FALSE)
      Grassberger <- bcShannon(Ns, Correction="Grassberger", CheckArguments=FALSE)
      return(max(ChaoShen, Grassberger))
    } else {
      if (Correction == "ZhangGrabchak") {
        # Weights
        w_v <- 1/V
        return(sum(Ps*sapply(1:length(Ns), S_v)))
      } else {
        return (bcShannon(Ns, Correction=Correction, CheckArguments=FALSE)) 
      }
    }
  }
  
  # Not Shannon
  if (Correction == "ZhangGrabchak") {
    # Weights
    i <- 1:(N-1)
    w_vi <- (i-q)/i
    w_v <- cumprod(w_vi)
    return(sum(Ps*sapply(1:length(Ns), S_v))/(1-q))
  }
  SampleCoverage <- Coverage(Ns, CheckArguments=FALSE)
  if (Correction == "ChaoShen" | Correction == "Best" | Correction == "HTGrassberger") {
    CPs <- SampleCoverage*Ns/N 
  } 
  if (Correction == "ChaoShen" | Correction == "Best") {
  	ChaoShen <- -sum(CPs^q * lnq(CPs,q) /(1 - (1-CPs)^N))
  } 
  if (Correction == "ChaoShen") {
    return(ChaoShen)
  }
  if (Correction == "Grassberger" | Correction == "Best") {
    Grassberger <- (1-N^(-q)*sum(Enq(Ns, q)))/(q-1)
  }
  if (Correction == "Grassberger") {
    return(Grassberger)
  }
  if (Correction == "Best") {
    return(max(ChaoShen, Grassberger))
  }
  if (Correction == "Holste") {
   return(1/(1-q)*(beta(length(Ns)+N, q)*sum(1/beta(Ns+1, q))-1))  
  } 
  if (Correction == "Bonachela") {
  return(1/(1-q)*(beta(2+N, q)*sum(1/beta(Ns+1, q))-1))  
  } 
  warning("Correction was not recognized")
  return (NA)
}
