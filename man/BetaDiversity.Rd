\name{BetaDiversity}
\alias{BetaDiversity}
\title{
  Unbiased beta diversity of a metacommunity
}
\description{
  Calculates the unbiased beta diversity of order \eqn{q} between communities.
}
\usage{
BetaDiversity(MC, q, Correction = "Best", Tree = NULL, Normalize = TRUE, 
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
    If \code{TRUE} (default), diversity is not affected by the height of the tree.\cr
    If \code{FALSE}, diversity is proportional to the height of the tree.
  }
  \item{CheckArguments}{
    Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  Entropy is calculated by \code{\link{BetaEntropy}} and transformed into diversity.
  
  Diversity values of communities are not defined: community entropies are averaged to obtain the metacommunity entropy wich is transformed into diversity (Marcon et al., in revision).
}
\value{
  An \code{\link{MCdiversity}} object containing diversity value of the metacommunity.
}
\references{
  Marcon, E., Scotti, I., Herault, B., Rossi, V. and Lang, G. (in revision). Generalization of the partitioning of Shannon diversity. \emph{PLOS One}.
  
  Marcon, E., Herault, B. (in prep). Partitioning phylodiversity. \emph{Methods in Ecology and Evolution}.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{BetaEntropy}}
}
\examples{
  # Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
  data(Paracou618)
  # Estimate Shannon beta diversity
  summary(BetaDiversity(Paracou618.MC, 1))
  # Compare without correction
  summary(BetaDiversity(Paracou618.MC, 1, Correction = "None"))
  # Estimate phylogenetic Shannon beta diversity
  summary(BetaDiversity(Paracou618.MC, 1, Tree = Paracou618.Taxonomy) -> e)
}
