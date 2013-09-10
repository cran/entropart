bcPhyloDiversity <-
function(Ns, q, Tree, Normalize = TRUE, Correction = "Best", CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  ppTree <- Preprocess.Tree(Tree)
  if (Normalize) {
    Height <- 1
  } else {
    Height <- ppTree$Height
  }  

  # Calculate normalized entropy (Height will be considered later)
  Diversity <- bcPhyloEntropy(Ns, q, ppTree, Normalize=TRUE, Correction, CheckArguments=FALSE)
  # Transform it into diversity
  Diversity$Cuts <- expq(Diversity$Cuts, q)
  Diversity$Total <- expq(Diversity$Total, q) * Height
  # Complete it
  Diversity$Function <- "bcPhyloDiversity" 
  Diversity$Distribution <- deparse(substitute(Ps))
  Diversity$Tree <- deparse(substitute(Tree))
  Diversity$Type <- "alpha or gamma"
  Diversity$Order <- q
  Diversity$Correction <- Correction
  
  class(Diversity) <- c("PhyloDiversity", "PhyloValue")
  
  return(Diversity)  
}