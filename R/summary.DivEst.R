summary.DivEst <-
function(object, ...) {
  
  cat("Diversity partitioning of order", object$Order, "of MetaCommunity", object$MetaCommunity)
  if (!object$Biased)  
    cat(" with correction:", object$Correction)
  cat("\n")
  
  if (!is.null(object$Tree)) {
    cat("Phylogenetic or functional diversity was calculated according to the tree", object$Tree, "\n")
    cat("Diversity is", ifelse(object$Normalized, "", "not"), "normalized\n")
  }
  cat("Alpha diversity of communities:", "\n")
  print(object$CommunityAlphaDiversities)
  cat("Total alpha diversity of the communities:", "\n")
  print(object$TotalAlphaDiversity)
  cat("Beta diversity of the communities:", "\n")
  print(object$TotalBetaDiversity)
  cat("Gamma diversity of the metacommunity:", "\n")
  print(object$GammaDiversity)
  
  cat("Quantiles of simulations (alpha, beta and gamma diviersity):\n")
  quant <- c(0, 0.1, 0.5, 0.1, seq(0.25, 0.75, 0.25), 0.9, 0.95, 0.99, 1)
  print(quantile(object$SimulatedDiversity["Alpha", ], probs = quant))
  print(quantile(object$SimulatedDiversity["Beta", ], probs = quant))
  print(quantile(object$SimulatedDiversity["Gamma", ], probs = quant))
  
  return(invisible(NULL))
}