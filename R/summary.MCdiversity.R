summary.MCdiversity <-
function(object, ...) {
  
  cat(object$Type, "diversity of order", object$Order, "of MetaCommunity", object$MetaCommunity, "with correction:", object$Correction, "\n")
  if (!is.null(object$Tree)) {
    cat("Phylogenetic or functional diversity was calculated according to the tree", object$Tree, "\n")
    cat("Diversity is", ifelse(object$Normalized, "", "not"), "normalized\n")
  }
  if (!is.null(object$Communities)) {
    cat("Diversity of communities:", "\n")
    print(object$Communities)
  }
  cat("Average diversity of the communities:", "\n")
  print(object$Total)

  return(invisible(NULL))
}