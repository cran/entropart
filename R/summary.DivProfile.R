summary.DivProfile <-
function(object, ...) {
  
  cat("Diversity profile of MetaCommunity", object$MetaCommunity)
  if (!object$Biased)  
    cat(" with correction:", object$Correction)
  cat("\n")
  
  if (!is.null(object$Tree)) {
    cat("Phylogenetic or functional diversity was calculated according to the tree", object$Tree, "\n")
    cat("Diversity is", ifelse(object$Normalized, "", "not"), "normalized\n")
  }
  
  cat("Diversity against its order:\n")
  Values <- cbind(object$Order, object$TotalAlphaDiversity, object$TotalBetaDiversity, object$GammaDiversity)
  colnames(Values) <- c("Order", "Alpha Diversity", "Beta Diversity", "Gamma Diversity")
  print(Values)
    
  return(invisible(NULL))
}