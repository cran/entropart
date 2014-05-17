summary.MCentropy <-
function(object, ...) {
  
  cat(object$Method, object$Type, "entropy of order", object$Order, "of metaCommunity", object$MetaCommunity, "with correction:", object$Correction, "\n", fill=TRUE)
  if (!is.null(object$Tree)) {
    cat("Phylogenetic or functional entropy was calculated according to the tree", object$Tree, "\n", fill=TRUE)
    cat("Entropy is", ifelse(object$Normalized, "", "not"), "normalized", fill=TRUE)
  }
  if (!is.null(object$Z)) {
    cat("Phylogenetic or functional entropy was calculated according to the similarity matrix", object$Z, "\n", fill=TRUE)
  }
  cat("Entropy of communities:", "\n")
  print(object$Communities)
  cat("Average entropy of the communities:", "\n")
  print(object$Total)

  return(invisible(NULL))
}