BetaDiversity <-
function(MC, q, Correction = "Best", Tree = NULL, Normalize = TRUE, CheckArguments = TRUE) 
{
  if (CheckArguments) 
    CheckentropartArguments()
  
  ppTree <- Preprocess.Tree(Tree)
  BetaEntropy  <- BetaEntropy (MC=MC, q=q, Correction=Correction, Tree=ppTree, Normalize=Normalize, CheckArguments=FALSE)
  AlphaEntropy <- AlphaEntropy(MC=MC, q=q, Correction=Correction, Tree=ppTree, Normalize=Normalize, CheckArguments=FALSE)

  Diversity <- list(
    MetaCommunity = deparse(substitute(MC)),
    Type = "beta",
    Order = q,
    Correction = Correction,
    Normalized = Normalize,
    Weights = MC$Wi, 
    Total = expq(BetaEntropy$Total / (1 - (q-1)*AlphaEntropy$Total), q)
  )
  if(!is.null(Tree))
    Diversity$Tree <- deparse(substitute(Tree)) 
  class(Diversity) <- "MCdiversity"
  
  return(Diversity)  
  
}
