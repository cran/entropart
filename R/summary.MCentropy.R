summary.MCentropy <-
function(object, ...) {
  
  cat(object$Type, "Entropy of order", object$Order, "of MetaCommunity", object$MetaCommunity, "with correction:", object$Correction, "\n")
  if (!is.null(object$Tree)) {
    cat("Phylogenetic or functional entropy was calculated according to the tree $Tree", "\n")
    cat("Entropy is", ifelse(object$Normalized, "", "not"), "normalized")
  }
  cat("Entropy of communities:", "\n")
  print(object$Communities)
  cat("Average entropy of the communities:", "\n")
  print(object$Total)

  return(invisible(NULL))
}