Originality.Species <-
function(Tree, q = 2, CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  # Tree must be either a phylog or a hclust object
  if (inherits(Tree, "phylog")) {
    phyTree <- Tree
  } else {
    if (inherits(Tree, "hclust")) {
      # build a phylog object
      phyTree <- ade4::hclust2phylog(Tree)
    } else {
      stop("Tree must be an object of class phylog or hclust")
    }
  }
 
  # Calculate the originality according to Rao entropy (exact computation)
  Originality <- ade4::originality(phyTree)
  Freq <- as.vector(t(Originality))
  names(Freq) <- rownames(Originality)

  
  # For q!=2, numerical optimization is required. Freq is used as a starting point.
  if (q != 2) {
    # bcPhyloEntropy accepts real values
    Optim <- stats::optim(Freq, 
                   function(Ns, q, Tree) bcPhyloEntropy(Ns, q, Tree, Correction="None", Normalize=FALSE, CheckArguments=FALSE)$Total, 
                   gr=NULL, q=q, Tree=Tree, lower=rep(0, length(Freq)), 
                   upper=rep(1, length(Freq)), method="L-BFGS-B", control=list(fnscale=-1))
    # Normalize
    Freq <-  Optim$par/sum(Optim$par)
  }
  return (Freq)
}
