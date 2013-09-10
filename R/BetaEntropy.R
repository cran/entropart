BetaEntropy <-
function(MC, q, Correction = "Best", Tree = NULL, Normalize = TRUE, CheckArguments = TRUE) 
{
  if (CheckArguments) 
    CheckentropartArguments()
    
  # Communities
  if (is.null(Tree)) {
    Communities <- apply(MC$Nsi, 2, bcTsallisBeta, Nexp=MC$Ns, q=q, Correction=Correction, CheckArguments=FALSE)
  } else {
    # Get the entropy of communities. Unlist to get a vector
    DetailedCommunities <- apply(MC$Nsi, 2, bcPhyloBetaEntropy, Nexp=MC$Ns, q=q, Tree=Tree, Normalize=Normalize, Correction=Correction, CheckArguments=FALSE)
    # Get $Total in each community's list
    Communities <- unlist(lapply(DetailedCommunities, function(x) x$Total))
  }
  # Weighted sum of Communities
  Total <- sum(Communities*MC$Wi)
  
  Entropy <- list(
    MetaCommunity = deparse(substitute(MC)),
    Type = "beta",
    Order = q,
    Correction = Correction,
    Normalized = Normalize,
    Weights = MC$Wi, 
    Communities = Communities, 
    Total=Total
  )
  if(!is.null(Tree))
    Entropy$Tree <- deparse(substitute(Tree)) 
  class(Entropy) <- "MCentropy"
  
  return(Entropy)  
}
