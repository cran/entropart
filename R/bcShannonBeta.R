bcShannonBeta <-
function(Ns, Nexp, Correction = "Best", CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  if (Correction == "ZhangGrabchak") {
    # Eliminate 0
    Y <- Nexp[Ns > 0]
    X <- Ns[Ns > 0]
    m <- sum(Nexp)
    n <- sum(Ns)
    Ps <- X/n
    V1 <- 1:m
    V2 <- 1:(n-1)
    # p_V_Ns1 is an array, containing (1 - Y_s/(m-j+1)) for each species (lines) and all j from 1 to m (Y_s may be 0)
    p_V_Ns1 <- outer(Y, V1, function(Y, j) 1- Y/(m-j+1))
    # p_V_Ns2 is an array, containing (1 - (X_s-1)/(n-j)) for each species (lines) and all j from 1 to n-1
    p_V_Ns2 <- outer(X, V2, function(X, j) 1- (X-1)/(n-j))
    # Useful values are products from j=1 to v, so prepare cumulative products
    p_V_Ns1 <- t(apply(p_V_Ns1, 1, cumprod))
    p_V_Ns2 <- t(apply(p_V_Ns2, 1, cumprod))
    # Sum of products weighted by w_v=1/v
    w_v1 <- 1/V1
    w_v2 <- 1/V2
    S_v <- function(s) {
      Usedv1 <- 1:(m-Y[s])
      Usedv2 <- 1:(n-X[s])
      return(sum(w_v1[Usedv1]*p_V_Ns1[s, Usedv1]) - sum(w_v2[Usedv2]*p_V_Ns2[s, Usedv2]))
    }
    return(sum(Ps*sapply(1:length(Ps), S_v)))
  } else {
    return (bcTsallisBeta(Ns, Nexp, 1, Correction, CheckArguments=FALSE))    
  }
}
