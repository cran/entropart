PhyloDiversity <-
function(Ps, q = 1, Tree, Normalize = TRUE, CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  ppTree <- Preprocess.Tree(Tree)
  if (Normalize) {
    Height <- 1
  } else {
    Height <- ppTree$Height
  }  
  # Calculate entropy
  Diversity <- PhyloEntropy(Ps, q, ppTree, Normalize=TRUE, CheckArguments=FALSE)
  # Transform it into diversity
  Diversity$Cuts <- expq(Diversity$Cuts, q)
  Diversity$Total <- expq(Diversity$Total, q) * Height
  # Complete it
  Diversity$Function <- "PhyloDiversity" 
  Diversity$Distribution <- deparse(substitute(Ps))
  Diversity$Tree <- deparse(substitute(Tree))
  Diversity$Type <- "alpha or gamma"
  Diversity$Order <- q
  
  class(Diversity) <- c("PhyloDiversity", "PhyloValue")

  return(Diversity)  
}
