summary.PhyloEntropy <-
function(object, ...) {
  
  cat(object$Type, "phylogenetic or functional entropy of order", object$Order, "of distribution", object$Distribution, fill=TRUE)
  if (!is.null(object$Correction)) {
    cat(" with correction:", object$Correction)
  }
  if (!is.null(object$Tree)) {
    cat("\nPhylogenetic or functional entropy was calculated according to the tree", object$Tree, "\n", fill=TRUE)
    cat("Entropy is", ifelse(object$Normalized, "", "not"), "normalized", fill=TRUE)
  }
  cat("\nEntropy equals:", object$Total)
  return(invisible(NULL))
}