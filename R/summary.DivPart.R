summary.DivPart <-
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
  
  return(invisible(NULL))
}