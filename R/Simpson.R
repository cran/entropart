Simpson <-
function(NorP, Correction="Lande", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  UseMethod("Simpson")
}


Simpson.ProbaVector <-
function(NorP, Correction="Lande", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  return (Tsallis(NorP, q=2, CheckArguments=FALSE))
}


Simpson.AbdVector <-
function(NorP, Correction="Lande", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  return(bcSimpson(Ns=NorP, Correction=Correction, CheckArguments=CheckArguments))
}


Simpson.integer <-
function(NorP, Correction="Lande", CheckArguments = TRUE, Ps = NULL, Ns = NULL)
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP or Ns must be provided.")
    }
  }
  return(bcSimpson(Ns=NorP, Correction=Correction, CheckArguments=CheckArguments))
}


Simpson.numeric <-
function(NorP, Correction="Lande", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
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
  
  if (abs(sum(NorP) - 1) < 3*.Machine$double.eps) {
    # Probabilities sum to 1, allowing rounding error
    return(Simpson.ProbaVector(NorP, CheckArguments=CheckArguments))
  } else {
    # Abundances
    return(Simpson.AbdVector(NorP, Correction=Correction, CheckArguments=CheckArguments))
  }
}


bcSimpson <-
function(Ns, Correction="Lande", CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  if (Correction == "Lande" | Correction == "Best") {
    N <- sum(Ns)
    return(N/(N-1)*Tsallis(as.ProbaVector(Ns), 2, CheckArguments=FALSE))  
  } else {
    return(bcTsallis(Ns, q=2, Correction, CheckArguments=FALSE)) 
  }
}
