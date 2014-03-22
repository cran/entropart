summary.MetaCommunity <-
function(object, ...) {
  
  cat("Meta-community (class 'MetaCommunity') made of ", object$N, "individuals in", object$Ncommunities, "communities and", object$Nspecies, "species.", "\n", fill=TRUE)
  cat(paste("Its sample coverage is ", object$SampleCoverage, "\n"), fill=TRUE)
  cat("Community weights are:", "\n")
  print(object$Wi)
  cat("Community sample numbers of individuals are:", "\n")
  print(object$Ni)
  cat("Community sample coverages are:", "\n")
  print(object$SampleCoverage.communities)

  return(invisible(NULL))
}