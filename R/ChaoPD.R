ChaoPD <-
function(Ps, q = 1, PhyloTree, Normalize = TRUE, CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  # Tree must be either a phylog, phylo or a hclust object
  if (inherits(PhyloTree, "phylog")) {
    # build a phylo object. Go through an hclust because as.phylo.phylog generates errors
    hTree <- stats::hclust(PhyloTree$Wdist^2/2, "average")
    Tree <- ape::as.phylo.hclust(hTree)
  } else {
    if (inherits(PhyloTree, "hclust")) {
      # build a phylo object
      Tree <- ape::as.phylo.hclust(PhyloTree)
    } else {
      if (inherits(PhyloTree, "phylo")) {
        Tree <- PhyloTree
      } else {
        if (inherits(PhyloTree, "PPtree")) {
          Tree <- PhyloTree$phyTree
        } else {
          stop("Tree must be an object of class phylo, phylog, hclust or PPtree")
        }
      }
    }
  }
  
  # Ps must be named
  if (is.null(names(Ps)))
    stop("Ps must be named to match the tree's species")
  # Verifiy all species are in the tree
  SpeciesNotFound <- setdiff(names(Ps), Tree$tip.label)
  if (length(SpeciesNotFound) > 0) {
    stop(paste("Species not found in the tree: ", SpeciesNotFound, collapse = "; "))
    print(SpeciesNotFound)
  }

  # Branch lengths
  Lengths <- Tree$edge.length
  # Get unnormalized probabilities p(b)
  ltips <- sapply(Tree$edge[, 2], function(node) geiger::tips(Tree, node))
  Branches <- unlist(lapply(ltips, function(VectorOfTips) sum(Ps[VectorOfTips])))
  # Normalize to get l(b)
  Tbar <- sum(Lengths*Branches)
  Lengths <- Lengths/Tbar
  
  # Eliminate zeros
  Lengths <- Lengths[Branches != 0]
  Branches <- Branches[Branches != 0]

  # Return normalized diversity
  if (q == 1) {
    return (prod(Branches^(-Lengths*Branches))*ifelse(Normalize, 1, Tbar))
  } else {
    return ((sum(Lengths*Branches^q))^(1/(1-q))*ifelse(Normalize, 1, Tbar))
  }
    
}
