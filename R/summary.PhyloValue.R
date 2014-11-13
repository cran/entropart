summary.PhyloValue <-
function(object, ...) {
  
  cat(object$Function, "applied to", object$Distribution, "along the tree", object$Tree, fill=TRUE)
  cat("\nResults are", ifelse(object$Normalized, "", "not"), "normalized", fill=TRUE)
  cat("\nThe average value is:", object$Total)
  cat("\n\nValues along the tree are:\n")
  print(object$Cuts)

  return(invisible(NULL))
}