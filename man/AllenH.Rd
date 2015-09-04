\name{AllenH}
\alias{AllenH}
\title{
  Phylogenetic Entropy of a Community
}
\description{
  Calculates the phylogenetic diversity of order \eqn{q} of a probability vector.
}
\usage{
AllenH(Ps, q = 1, PhyloTree, Normalize = TRUE, CheckArguments = TRUE)
}
\arguments{
  \item{Ps}{
    A probability vector, summing to 1.
  }
  \item{q}{
    A number: the order of entropy. Default is 1.
  }
  \item{PhyloTree}{
    An object of class \code{\link{hclust}} or \code{\link{phylog}}. The tree is not necessarily ultrametric. Computation is must faster if the tree is \code{\link{phylog}} or \code{\link{PPtree}}
  }
  \item{Normalize}{
    If \code{TRUE} (default), diversity is not affected by the height of the tree.\cr
    If \code{FALSE}, it is proportional to the height of the tree.\cr
  }
  \item{CheckArguments}{
    Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  The phylogenetic entropy is calculated following Allen \emph{et al.} (2009) for order \eqn{q=1} and Leinster and Cobold (2011) for other orders.The result is identical to the total entropy calculated by \code{\link{PhyloEntropy}} but it is much faster. A single value is returned instead of a \code{\link{PhyloEntropy}} object, and no bias correction is available.
  
  The \code{Normalize} argument allows normalizing entropy by the height of the tree, similarly to \code{\link{ChaoPD}}.
  
  This function allows very efficient calculation when the \code{Tree} is of class \code{\link{phylog}} or \code{\link{PPtree}} (\code{\link{hclust}} trees must be converted to \code{\link{phylog}}, a slow process).
  
  Diversity can be calculated for non ultrametric trees following Leinster and Cobold (2011) even though the meaning of the result is not so clear.
}
\value{
  A number equal the entropy of the community.
}
\references{
  Allen, B., Kon, M. and Bar-Yam, Y. (2009). A New Phylogenetic Diversity Measure Generalizing the Shannon Index and Its Application to Phyllostomid Bats. \emph{American Naturalist} 174(2): 236-243.
  
  Leinster, T. and Cobbold, C. (2011). Measuring diversity: the importance of species similarity. \emph{Ecology} 93(3): 477-489.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{PhyloEntropy}}, \code{\link{ChaoPD}}
}
\examples{
  # Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest
  #      and their taxonomy)
  data(Paracou618)
  # Ps is the vector of probabilities
  Ps <- as.ProbaVector(Paracou618.MC$Ns)
  
  # Calculate the phylogenetic Shannon diversity of the plot
  AllenH(Ps, 1, Paracou618.Taxonomy, Normalize=TRUE)
  
  # Calculate it using PhyloEntropy: more powerful but much slower is the tree has many periods
  PhyloEntropy(Ps, 1, Paracou618.Taxonomy, Normalize=TRUE) -> phyE
  summary(phyE)
}
