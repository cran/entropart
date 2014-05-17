bcPhyloEntropy <-
function(Ns, q = 1, Tree, Normalize = TRUE, Correction = "Best", CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  # Calculate the PhyloValue
  Entropy <- PhyloApply(Tree, bcTsallis, Ns, Normalize, q=q, Correction=Correction, CheckArguments=FALSE)
  # Complete it
  Entropy$Function <- "PhyloEntropy" 
  Entropy$Distribution <- deparse(substitute(Ns)) 
  Entropy$Tree <- deparse(substitute(Tree))
  Entropy$Type <- "alpha or gamma"
  Entropy$Order <- q
  Entropy$Correction <- Correction
  
  class(Entropy) <- c("PhyloEntropy", class(Entropy))
  
  return (Entropy)                         
}
