\name{ChaoPD}
\alias{ChaoPD}
\title{
  Phylogenetic Diversity of a Community
}
\description{
  Calculates the phylogenetic diversity of order \eqn{q} of a probability vector.
}
\usage{
ChaoPD(Ps, q = 1, PhyloTree, Normalize = TRUE, CheckArguments = TRUE)
}
\arguments{
  \item{Ps}{
    A probability vector, summing to 1.
  }
  \item{q}{
    A number: the order of diversity. Default is 1.
  }
  \item{PhyloTree}{
    An object of class \code{\link{hclust}}, \code{\link{phylo}}, \code{\link{phylog}} or \code{\link{PPtree}}. The tree is not necessarily ultrametric.
  }
  \item{Normalize}{
    If \code{TRUE} (default), diversity is not affected by the height of the tree.\cr
    If \code{FALSE}, it is proportional to the height of the tree.
  }
  \item{CheckArguments}{
    Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  The phylogenetic diversity is calculated following Chao \emph{et al.} (2010). The result is identical to the total diversity calculated by \code{\link{PhyloDiversity}} but it is much faster. A single value is returned instead of a \code{\link{PhyloDiversity}} object, and no bias correction is available.
  
  The \code{Normalize} arguments allows calculating either \eqn{^{q}\bar{D}(T)} (if \code{TRUE}) or \eqn{^{q}PD(T)} if \code{FALSE}.

  Diversity can be calculated for non ultrametric trees following Chao \emph{et al.} (2010) even though the meaning of the result is not so clear (Leinster and Cobold, 2011).
}
\value{
  A number equal the diversity of the community.
}
\references{
  Chao, A., Chiu, C.-H. and Jost, L. (2010). Phylogenetic diversity measures based on Hill numbers. \emph{Philosophical Transactions of the Royal Society B} 365(1558): 3599-609.
  
  Leinster, T. and Cobbold, C. (2011). Measuring diversity: the importance of species similarity. \emph{Ecology} 93(3): 477-489.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{PhyloDiversity}}, \code{\link{AllenH}}
}
\examples{
  # Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest
  #      and their taxonomy)
  data(Paracou618)
  # Ps is the vector of probabilities
  Ps <- Paracou618.MC$Ps
  
  # Calculate the phylogenetic Simpson diversity of the plot
  ChaoPD(Paracou618.MC$Ps, 2, Paracou618.Taxonomy, Normalize=TRUE) -> ChaoD
  ChaoD
  
  # Calculate it using PhyloDiversity
  # (more powerful but much slower is the tree has many periods)
  PhyloDiversity(Paracou618.MC$Ps, 2, Paracou618.Taxonomy, Normalize=TRUE) -> phyD
  summary(phyD)
}
