Preprocess.Tree <-
function(Tree)
{
  # The tree may be NULL or already processed
  if (is.null(Tree) | inherits(Tree, "PPtree")) return (Tree)

  # Tree must be either a phylog or a hclust object
  if (inherits(Tree, "phylog")) {
    phyTree <- Tree
    # build an hclust object to use cutree later. Distances in $Wdist are actually 2*sqrt(distance)
    hTree <- hclust(Tree$Wdist^2/2, "average")
  } else {
    if (inherits(Tree, "hclust")) {
      hTree <- Tree
      # build a phylog object to use $droot later
      phyTree <- ade4::hclust2phylog(Tree)
    } else {
      stop("Tree must be an object of class phylog or hclust")
    }
  }

  # Calculate distances between nodes and leaves ($droot are distances from root)
  DistancesFromLeaves <- max(phyTree$droot)-phyTree$droot
  # Get a sorted list of cuts (eliminate leaves)
  Cuts <- sort(DistancesFromLeaves[setdiff(names(DistancesFromLeaves), names(phyTree$leaves))])
  # Calculate intervals between cuts (add 0 to Cuts to get the first interval)
  Intervals <- diff(c(0, Cuts))
  # Eliminate 0 intervals (when a node contains more than 2 tips), including rounding errors
  RoundingError <- max(phyTree$droot) * 10 * .Machine$double.eps
  Cuts <- Cuts[Intervals > RoundingError]
  Intervals <- Intervals[Intervals > RoundingError]

  PPtree <- list(
    phyTree   = phyTree,
    hTree     = hTree,
    Height    = Cuts[length(Cuts)],
    Cuts      = Cuts,
    Intervals = Intervals
    )
  class(PPtree) <- "PPtree"
  return (PPtree)
}
