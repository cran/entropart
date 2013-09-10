bcPhyloBetaEntropy <-
function(Ns, Nexp, q, Tree, Normalize = TRUE, Correction = "Best", CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()

  # Prepare NorP
  NandNexp <- matrix(c(Ns, Nexp), nrow = length(Ns), ncol = 2, dimnames = list(names(Ns), c("Ns", "Nexp")))
  #  Calculate the PhyloValue. Intermediate function is necessary to separate N and Nexp before calling bcTsallisBeta
  Entropy <- PhyloApply(Tree, function(NandNexp, q, Correction, CheckArguments) bcTsallisBeta(NandNexp[, "Ns"], NandNexp[, "Nexp"], q, Correction, CheckArguments), NandNexp, Normalize, q=q, Correction=Correction, CheckArguments=FALSE)
  # Complete it
  Entropy$Function <- "PhyloEntropy" 
  Entropy$Distribution <- deparse(substitute(Ps)) 
  Entropy$Tree <- deparse(substitute(Tree))
  Entropy$Type <- "beta"
  Entropy$Order <- q
  Entropy$Correction <- Correction
  
  class(Entropy) <- c("PhyloEntropy", class(Entropy))
  
  return (Entropy)
}
