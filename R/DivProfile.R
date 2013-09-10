DivProfile <-
function(q.seq = seq(0, 2, .1), MC, Biased = TRUE, Correction = "Best", Tree = NULL, Normalize = TRUE, CheckArguments = TRUE) 
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

  # Calculate diversity profile
  Diversity.seq <- sapply(q.seq, DivPart, MC = MC, Biased = Biased, Correction = Correction, Tree = ppTree, Normalize = Normalize, CheckArguments =FALSE)
  # Prepare a matrix of results
  DivProfile <- list(MetaCommunity = unlist(Diversity.seq["MetaCommunity", 1], use.names=FALSE),
                     Order = unlist(Diversity.seq["Order", ], use.names=FALSE), 
                     Biased = unlist(Diversity.seq["Biased", 1], use.names=FALSE), 
                     Correction = unlist(Diversity.seq["Correction", 1], use.names=FALSE),
                     Normalized = unlist(Diversity.seq["Normalized", 1], use.names=FALSE),
                     TotalAlphaDiversity = unlist(Diversity.seq["TotalAlphaDiversity", ]), 
                     TotalBetaDiversity = unlist(Diversity.seq["TotalBetaDiversity", ]), 
                     GammaDiversity = unlist(Diversity.seq["GammaDiversity", ]), 
                     TotalAlphaEntropy =  unlist(Diversity.seq["TotalAlphaEntropy", ]), 
                     TotalBetaEntropy = unlist(Diversity.seq["TotalBetaEntropy", ]), 
                     GammaEntropy =  unlist(Diversity.seq["GammaEntropy", ])
                    )
  if(!is.null(Tree))
    DivProfile$Tree <- deparse(substitute(Tree)) 
  class(DivProfile) <- "DivProfile"
    
  return (DivProfile)
}
