\name{BetaEntropy}
\alias{BetaEntropy}
\title{
  Unbiased beta entropy of a metacommunity
}
\description{
  Calculates the unbiased beta entropy of order \eqn{q} between communities.
}
\usage{
BetaEntropy(MC, q, Correction = "Best", Tree = NULL, Normalize = TRUE, 
  CheckArguments = TRUE)
}
\arguments{
  \item{MC}{
    A \code{\link{MetaCommunity}} object.
  }
  \item{q}{
    A number: the order of diversity.
  }
  \item{Correction}{
    A string containing one of the possible corrections accepted by \code{\link{bcTsallisBeta}}. \code{"Best"} is the default value.
  }
  \item{Tree}{
    An object of class \code{\link{hclust}} or \code{\link{phylog}}. The tree must be ultrametric.
  }
  \item{Normalize}{
    If \code{TRUE} (default), the entropy returned by the function is normalized by the height of the tree (it is the weighted average value of the entropy in each slice).\cr
    If \code{FALSE}, it is the unnormalized weighted sum of the results.
  }
  \item{CheckArguments}{
    Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  If \code{Tree} is \code{NULL}, then generalized entropy is calculated by \code{\link{bcTsallisBeta}}, else phylogenetic entropy is calculated by \code{\link{bcPhyloBetaEntropy}}.
  
  The unbiased beta entropy of each community is calculated and summed according to community weights.

  Note that beta entropy is related to alpha entropy (if \eqn{q} is not 1) and cannot be compared accross communities (Jost, 2007). Do rather calculate the \code{\link{BetaDiversity}} of the metacommunity.
}
\value{
  An \code{\link{MCentropy}} object containing entropy values of each community and of the metacommunity.
}
\references{
  Marcon, E., Scotti, I., Herault, B., Rossi, V. and Lang, G. (in revision). Generalization of the partitioning of Shannon diversity. \emph{PLOS One}.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{bcTsallisBeta}}, \code{\link{BetaDiversity}}
}
\examples{
  # Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
  data(Paracou618)
  # Estimate Shannon beta entropy
  summary(BetaEntropy(Paracou618.MC, 1))
  # Compare without correction
  summary(BetaEntropy(Paracou618.MC, 1, Correction = "None"))
  # Estimate phylogenetic Shannon beta entropy
  summary(BetaEntropy(Paracou618.MC, 1, Tree = Paracou618.Taxonomy) -> e)
  plot(e)  
}
