HurlbertD <-
function(NorP, k = 2, CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  UseMethod("HurlbertD")
}


HurlbertD.ProbaVector <-
function(NorP, k = 2, CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  # Calculate Hurlbert's index
  Sk <- Hurlbert(NorP, k = k, CheckArguments = FALSE)
  # Find the effective number of species numerically
  f <- function(D, S, k) D*(1-(1-1/D)^k)-S
  Dk <- stats::uniroot(f, c(1, 1E+7), S=Sk, k=k)
  return (Dk$root)
}

HurlbertD.AbdVector <-
function(NorP, k = 2, CheckArguments = TRUE, Ps = NULL, Ns = NULL)
{
  return(bcHurlbertD(Ns=NorP, k=k, CheckArguments=CheckArguments))
}


HurlbertD.integer <-
function(NorP, k = 2, CheckArguments = TRUE, Ps = NULL, Ns = NULL)
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP Ns must be provided.")
    }
  }
  return(bcHurlbertD(Ns=NorP, k=k, CheckArguments=CheckArguments))
}


HurlbertD.numeric <-
function(NorP, k = 2, CheckArguments = TRUE, Ps = NULL, Ns = NULL)
{
  if (missing(NorP)){
    if (!missing(Ps)) {
      NorP <- Ps
    } else {
      if (!missing(Ns)) {
        NorP <- Ns
      } else {
        stop("An argument NorP or Ps or Ns must be provided.")
      }
    }
  }
  
  if (abs(sum(NorP) - 1) < 10*.Machine$double.eps) {
    # Probabilities sum to 1, allowing rounding error
    return(HurlbertD.ProbaVector(NorP, k=k, CheckArguments=CheckArguments))
  } else {
    # Abundances
    return(HurlbertD.AbdVector(NorP, k=k, CheckArguments=CheckArguments))
  }
}


bcHurlbertD <-
function(Ns, k = 2, CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  # Calculate Hurlbert's index
  Sk <- bcHurlbert(Ns, k = k, CheckArguments = FALSE)
  # Find the effective number of species numerically
  f <- function(D, S, k) D*(1-(1-1/D)^k)-S
  Dk <- stats::uniroot(f, c(1, 1E+7), S=Sk, k=k)
  return (Dk$root)
}
