GammaEntropy <-
function(MC, q, Correction = "Best", Tree = NULL, Normalize = TRUE, CheckArguments = TRUE) 
{
  if (CheckArguments) 
    CheckentropartArguments()
  
  if (is.null(Tree)) {
    Entropy <- bcTsallis(MC$Ns, q, Correction)
  } else {
    Entropy <- bcPhyloEntropy(MC$Ns, q, Tree, Normalize, Correction)$Total
  }
  
  return(Entropy)
}
