AlphaDiversity <-
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
  AlphaEntropy <- AlphaEntropy(MC, q, Correction, ppTree, Normalize=TRUE)
  Diversity <- list(
    MetaCommunity = deparse(substitute(MC)),
    Type = "alpha",
    Order = q,
    Correction = Correction,
    Normalized = Normalize,
    Weights = MC$Wi, 
    Communities = expq(AlphaEntropy$Communities, q) * Height,
    Total = expq(AlphaEntropy$Total, q)* Height
    )
  if(!is.null(Tree))
    Diversity$Tree <- deparse(substitute(Tree)) 
  class(Diversity) <- "MCdiversity"
  
  return(Diversity)  
}
