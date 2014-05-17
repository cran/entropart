summary.MCdiversity <-
function(object, ...) {
  
  cat(object$Method, object$Type, "diversity of order", object$Order, "of metaCommunity", object$MetaCommunity, "with correction:", object$Correction, "\n", fill=TRUE)
  if (!is.null(object$Tree)) {
    cat("Phylogenetic or functional diversity was calculated according to the tree", object$Tree, "\n", fill=TRUE)
    cat("Diversity is", ifelse(object$Normalized, "", "not"), "normalized\n", fill=TRUE)
  }
  if (!is.null(object$Z)) {
    cat("Phylogenetic or functional entropy was calculated according to the similarity matrix", object$Z, "\n", fill=TRUE)
  }
  if (!is.null(object$Communities)) {
    cat("Diversity of communities:", "\n")
    print(object$Communities)
  }
  cat("Average diversity of the communities:", "\n")
  print(object$Total)

  return(invisible(NULL))
}