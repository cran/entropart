summary.PhyloEntropy <-
function(object, ...) {
  
  cat(object$Type, "phylogenetic or functional entropy of order", object$Order, "of distribution", object$Distribution)
  if (!is.null(object$Correction)) {
    cat(" with correction:", object$Correction)
  }
  if (!is.null(object$Tree)) {
    cat("\nPhylogenetic or functional entropy was calculated according to the tree", object$Tree, "\n")
    cat("Entropy is", ifelse(object$Normalized, "", "not"), "normalized")
  }
  cat("\nEntropy equals:", object$Total)
  return(invisible(NULL))
}