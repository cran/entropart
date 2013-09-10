PhyloEntropy <-
function(Ps, q, Tree, Normalize = TRUE, CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  # Calculate the PhyloValue
  Entropy <- PhyloApply(Tree, Tsallis, Ps, Normalize, q=q, CheckArguments=FALSE)
  # Complete it
  Entropy$Function <- "PhyloEntropy" 
  Entropy$Distribution <- deparse(substitute(Ps)) 
  Entropy$Tree <- deparse(substitute(Tree))
  Entropy$Type <- "alpha or gamma"
  Entropy$Order <- q
  
  class(Entropy) <- c("PhyloEntropy", class(Entropy))

  return (Entropy)
}
