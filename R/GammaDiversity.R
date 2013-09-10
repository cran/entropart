GammaDiversity <-
function(MC, q, Correction = "Best", Tree = NULL, Normalize = TRUE, CheckArguments = TRUE) 
{
  if (CheckArguments) 
    CheckentropartArguments()
  
  ppTree <- Preprocess.Tree(Tree)
  if (Normalize) {
    Height <- 1
  } else {
    Height <- ppTree$Height
  }  
  GammaEntropy <- GammaEntropy(MC, q, Correction, Tree, Normalize=TRUE)
  
  return(expq(GammaEntropy, q) * Height)  
}
