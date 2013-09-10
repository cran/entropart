summary.PhyloDiversity <-
function(object, ...) {
  
  cat(object$Type, "phylogenetic or functional diversity of order", object$Order, "of distribution", object$Distribution)
  if (!is.null(object$Correction)) {
    cat(" with correction:", object$Correction)
  }
  if (!is.null(object$Tree)) {
    cat("\nPhylogenetic or functional diversity was calculated according to the tree", object$Tree, "\n")
    cat("Diversity is", ifelse(object$Normalized, "", "not"), "normalized")
  }
  cat("\nDiversity equals:", object$Total)
  return(invisible(NULL))
}