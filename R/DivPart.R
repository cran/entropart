DivPart <-
function(q = 1, MC, Biased = TRUE, Correction = "Best", Tree = NULL, Normalize = TRUE, Z = NULL, CheckArguments = TRUE) 
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
    AlphaEntropy <- AlphaEntropy(MC, q, "None", ppTree, Normalize, Z, CheckArguments=FALSE)
    GammaEntropy <- GammaEntropy(MC, q, "None", ppTree, Normalize, Z, CheckArguments=FALSE)
    BetaEntropy  <- BetaEntropy (MC, q, "None", ppTree, Normalize, Z, CheckArguments=FALSE)
  } else {
    AlphaEntropy <- AlphaEntropy(MC, q, Correction, ppTree, Normalize, Z, CheckArguments=FALSE)
    GammaEntropy <- GammaEntropy(MC, q, Correction, ppTree, Normalize, Z, CheckArguments=FALSE)
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
  if(is.null(Z)) {
    DivPart$Method <- "HCDT"
  } else {
    DivPart$Method <- "Similarity-based"
    DivPart$Z <- deparse(substitute(Z))  
  }
  class(DivPart) <- "DivPart"
  
  return (DivPart)
}
