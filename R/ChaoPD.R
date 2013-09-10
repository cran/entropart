ChaoPD <-
function(Ps, q, PhyloTree, Normalize = TRUE, CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  # PhyloTree must be either a phylog or a hclust object
  if (inherits(PhyloTree, "PPtree")) {
    phyTree <- PhyloTree$phyTree
  } else {
    if (inherits(PhyloTree, "phylog")) {
      phyTree <- PhyloTree
    } else {
      if (inherits(PhyloTree, "hclust")) {
        phyTree <- hclust2phylog(PhyloTree)
      } else {
        stop("PhyloTree must be an object of class phylog or hclust")
      }
    }
  }
  
  # Verifiy all species are inthe tree
  SpeciesNotFound <- setdiff(names(Ps), names(phyTree$leaves))
  if (length(SpeciesNotFound) > 0) {
    stop(paste("Species not found in the tree: ", SpeciesNotFound, collapse = "; "))
    print(SpeciesNotFound)
  }
  
  
  # Prepare the vector of branches
  Names <- c(names(phyTree$leaves), names(phyTree$nodes))
  Branches <- vector("numeric", length(Names))
  names(Branches) <- Names
  # Eliminate "Root"
  Branches <- Branches[-length(Branches)]
  # Get unnormalized probabilities p(b)
  Branches[names(Ps)]=Ps
  for (NodeName in names(phyTree$nodes[-length(phyTree$nodes)])) {
    Branches[NodeName] <- sum(Branches[phyTree$parts[[NodeName]]])
  }
  
  # Lengths
  Lengths <- phyTree$droot[-length(phyTree$droot)]
  for (PartName in names(Lengths)) {
    Lengths[phyTree$parts[[PartName]]] <- phyTree$droot[phyTree$parts[[PartName]]]-phyTree$droot[PartName]
  }
  # Normalize to get l(b)
  Tbar <- mean(phyTree$droot[names(phyTree$leaves)])
  Lengths <- Lengths/Tbar

  # Return normalized diversity
  if (q == 1) {
    return (prod(Branches^(-Lengths*Branches))*ifelse(Normalize, 1, Tbar))
  } else {
    return ((sum(Lengths*Branches^q))^(1/(1-q))*ifelse(Normalize, 1, Tbar))
  }
    
}
