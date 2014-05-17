\name{PhyloEntropy}
\alias{PhyloEntropy}
\alias{bcPhyloEntropy}
\alias{is.PhyloEntropy}
\alias{summary.PhyloEntropy}
\title{
  Phylogenetic Entropy of a community
}
\description{
  Calculates the phylogenetic entropy of order \eqn{q} of a probability vector.
}
\usage{
PhyloEntropy(Ps, q = 1, Tree, Normalize = TRUE, CheckArguments = TRUE)
bcPhyloEntropy(Ns, q = 1, Tree, Normalize = TRUE, Correction = "Best", 
  CheckArguments = TRUE)
is.PhyloEntropy(x)
\method{summary}{PhyloEntropy}(object, \dots)
}
\arguments{
  \item{Ps}{
    A probability vector, summing to 1.
  }
  \item{Ns}{
    A numeric vector containing species abundances.
  }
  \item{q}{
    A number: the order of entropy. Default is 1.
  }
  \item{Tree}{
    An object of class \code{\link{hclust}} or \code{\link{phylog}}. The tree must be ultrametric.
  }
  \item{Normalize}{
    If \code{TRUE} (default), the \code{Total} entropy returned by the function is normalized by the height of the tree (it is the weighted average value of the entropy in each slice).\cr
    If \code{FALSE}, it is the unnormalized weighted sum of the results.
  }
  \item{Correction}{
    A string containing one of the possible corrections: \code{"None"} (no correction), \code{"ChaoShen"}, \code{"Grassberger"}, \code{"Holste"}, \code{"Bonachela"} or \code{"Best"}, the default value.
  }
  \item{CheckArguments}{
    Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
  \item{x}{
    An object to be tested or plotted
  }
  \item{object}{
    A \code{PhyloDiversity} object to be summarized.
  }
  \item{\dots}{
    Additional arguments to be passed to the generic methods.
  }
}
\details{
  The phylogenetic entropy is its generalization of HCDT entropy to unequal species distances (Pavoine et al., 2009).
  
  Calculation relies on \code{\link{Tsallis}} and \code{\link{PhyloApply}}.
  
  Intervals separate two cuts in a tree: no node is found at heights contained in an interval.
  
  Bias correction requires the number of individuals to estimate sample \code{\link{Coverage}}. Use \code{bcPhyloEntropy} and choose the \code{Correction}.
}
\value{
  An object of class \code{PhyloEntropy} is a list:
  \item{Distribution}{The distribution used to calculate entropy}
  \item{Function}{The function used to calculate entropy}
  \item{Tree}{The functional or phylogenetic tree used to calculate entropy}
  \item{Normalized}{Logical. Indicates whether phyloentropy is normalized or proportional to the height of the tree.}
  \item{Type}{The type of entropy ("alpha", "beta" or "gamma").}
  \item{Order}{The order of entropy \eqn{q}.}
  \item{Cuts}{A named vector containing values of neutral entropy along the tree. Names are cut ends, \emph{i.e.} the ends of intervals (the first interval starts at 0 for leaves, the max value is the height of the tree).}
  \item{Total}{A value equal the total entropy multiplied by the tree height if \code{Normalize} is FALSE.}
  
  \code{is.PhyloEntropy} returns \code{TRUE} if the object is of class \code{PhyloEntropy}.
  
  \code{summary.PhyloEntropy} returns a summary of the object's value.

  \code{PhyloEntropy} objects can be plotted by \code{\link{plot.PhyloValue}} because \code{PhyloEntropy} objects are also of class \code{PhyloValue}.
}
\references{
  Marcon, E., Herault, B. (2014). Decomposing  Phylodiversity. \emph{HAL} hal-00946177(version 1).
  
  Pavoine, S., Love, M. S. and Bonsall, M. B. (2009). Hierarchical partitioning of evolutionary and ecological patterns in the organization of phylogenetically-structured species assemblages: Application to rockfish (genus: Sebastes) in the Southern California Bight. \emph{Ecology Letters} 12(9): 898-908.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{Tsallis}}, \code{\link{PhyloDiversity}}
}
\examples{
  # Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest
  #      and their taxonomy)
  data(Paracou618)
  # Ps is the vector of probabilities
  Ps <- Paracou618.MC$Ps
  # Calculate the phylogenetic Shannon entropy of the plot
  summary(PhyloEntropy(Ps, 1, Paracou618.Taxonomy) -> e)
  plot(e)
  
  # Ns is the vector of abundances of the metacommunity
  Ns <- Paracou618.MC$Ns
  # Calculate the phylogenetic Shannon entropy of the plot
  summary(bcPhyloEntropy(Ns, 1, Paracou618.Taxonomy, Correction = "Best") -> e)
  plot(e)

}
