\name{GammaDiversity}
\alias{GammaDiversity}
\title{
  Unbiased gamma diversity of a metacommunity
}
\description{
  Calculates the unbiased diversity of order \eqn{q} of a metacommunity.
}
\usage{
GammaDiversity(MC, q, Correction = "Best", Tree = NULL, Normalize = TRUE, 
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
    A string containing one of the possible corrections: \code{"ChaoShen"}, \code{"Grassberger"}, \code{"Holste"}, \code{"Bonachela"} or \code{"Best"}, the default value.
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
  Entropy is calculated by \code{\link{GammaEntropy}} and transformed into diversity.
}
\value{
  The metacommunity's gamma entropy.
}
\references{
  Marcon, E., Scotti, I., Herault, B., Rossi, V. and Lang, G. (in revision). Generalization of the partitioning of Shannon diversity. \emph{PLOS One}.
  
  Marcon, E., Herault, B. (in prep). Partitioning phylodiversity. \emph{Methods in Ecology and Evolution}.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{GammaEntropy}}
}
\examples{
  # Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
  data(Paracou618)
  # Calculate Simpson gamma diversity
  GammaDiversity(Paracou618.MC, 2)
  # Compare without correction
  GammaDiversity(Paracou618.MC, 2, Correction = "None")
  # Estimate phylogenetic Simpson gamma diversity
  GammaDiversity(Paracou618.MC, 2, Tree = Paracou618.Taxonomy)   
}
