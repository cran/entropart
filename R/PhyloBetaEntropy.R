PhyloBetaEntropy <-
function(Ps, Pexp, q = 1, Tree, Normalize = TRUE, CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  # Prepare NorP
  PandPexp <- matrix(c(Ps, Pexp), nrow = length(Ps), ncol = 2, dimnames = list(names(Ps), c("Ps", "Pexp")))
  # Calculate the PhyloValue. Intermediate function is necessary to separate P and Pexp before calling TsallisBeta.P
  Entropy <- PhyloApply(Tree, function(PandPexp, q, CheckArguments) TsallisBeta(PandPexp[, "Ps"], PandPexp[, "Pexp"], q, CheckArguments), PandPexp, Normalize, q=q, CheckArguments=FALSE)
  # Complete it
  Entropy$Function <- "PhyloBetaEntropy" 
  Entropy$Distribution <- c(deparse(substitute(Ps)), "compared to", deparse(substitute(Pexp)))
  Entropy$Tree <- deparse(substitute(Tree))
  Entropy$Type <- "beta"
  Entropy$Order <- q
  
  class(Entropy) <- c("PhyloEntropy", class(Entropy))
  
  return (Entropy)
}
