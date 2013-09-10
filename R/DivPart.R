DivPart <-
function(q, MC, Biased = TRUE, Correction = "Best", Tree = NULL, Normalize = TRUE, CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  # Preprocess the tree
  ppTree <- Preprocess.Tree(Tree)
  if (Normalize) {
    Height <- 1
  } else {
    Height <- ppTree$Height
  }  

  # Alpha and beta entropy of communities
  if (Biased) {
    AlphaEntropy <- AlphaEntropy(MC, q, "None", ppTree, Normalize)
    GammaEntropy <- GammaEntropy(MC, q, "None", ppTree, Normalize)
    BetaEntropy  <- BetaEntropy (MC, q, "None", ppTree, Normalize)
  } else {
    AlphaEntropy <- AlphaEntropy(MC, q, Correction, ppTree, Normalize)
    GammaEntropy <- GammaEntropy(MC, q, Correction, ppTree, Normalize)
    # beta is calculated as gamma-alpha to ensure continuity. Community beta entropy is not calculated.
    BetaEntropy  <- list(Communities = NA, Total = GammaEntropy - AlphaEntropy$Total)      
  }
  # Total Diversities
  AlphaDiversity <- expq(AlphaEntropy$Total / Height, q) * Height
  BetaDiversity  <- expq(BetaEntropy$Total / Height / (1 - (q-1)*AlphaEntropy$Total/Height), q) 
  GammaDiversity <- expq(GammaEntropy / Height, q) * Height
  
  DivPart <- (list(
    MetaCommunity = deparse(substitute(MC)),
    Order = q, 
    Biased = Biased, 
    Correction = Correction,
    Normalized = Normalize,
    TotalAlphaDiversity = AlphaDiversity, 
    TotalBetaDiversity = BetaDiversity, 
    GammaDiversity = GammaDiversity, 
    CommunityAlphaDiversities = expq(AlphaEntropy$Communities / Height, q) * Height, 
    TotalAlphaEntropy = AlphaEntropy$Total, 
    TotalBetaEntropy = BetaEntropy$Total, 
    GammaEntropy = GammaEntropy, 
    CommunityAlphaEntropies = AlphaEntropy$Communities, 
    CommunityBetaEntropies = BetaEntropy$Communities
    ))
  if(!is.null(Tree))
    DivPart$Tree <- deparse(substitute(Tree)) 
  class(DivPart) <- "DivPart"
  
  return (DivPart)
}
