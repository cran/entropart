Rao <-
function(NorP, Tree, Correction = "Lande", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  UseMethod("Rao")
}


Rao.ProbaVector <-
function(NorP, Tree, Correction = "Lande", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  return(AllenH(NorP, q=2, PhyloTree=Tree, Normalize=FALSE, CheckArguments=FALSE))
}


Rao.AbdVector <-
function(NorP, Tree, Correction = "Lande", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  return(bcRao(Ns=NorP, Tree=Tree, Correction=Correction, CheckArguments=CheckArguments))
}


Rao.integer <-
function(NorP, Tree, Correction = "Lande", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP or Ns must be provided.")
    }
  }
  return(bcRao(Ns=NorP, Tree=Tree, Correction=Correction, CheckArguments=CheckArguments))
}


Rao.numeric <-
function(NorP, Tree, Correction = "Lande", CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
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
    return(Rao.ProbaVector(NorP, Tree=Tree, CheckArguments=CheckArguments))
  } else {
    # Abundances
    return(Rao.AbdVector(NorP, Tree=Tree, Correction=Correction, CheckArguments=CheckArguments))
  }
}


bcRao <-
function(Ns, Tree, Correction="Lande", CheckArguments = TRUE)
{
  if (CheckArguments)
    CheckentropartArguments()
  
  if (Correction == "Lande"  | Correction == "Best") {
    Nrecords <- sum(Ns)
    return(Nrecords/(Nrecords-1)*AllenH(Ns/sum(Ns), q=2, PhyloTree=Tree, Normalize=FALSE, CheckArguments=FALSE))  
  } else {
    return(bcPhyloEntropy(Ns, q = 2, Tree = Tree, Normalize = FALSE, Correction = Correction, CheckArguments = FALSE)) 
  }
}
