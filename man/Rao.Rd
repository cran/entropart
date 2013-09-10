\name{Rao}
\alias{Rao}
\title{
  Rao Quadratic Entropy of a Community
}
\description{
  Calculates Rao's quadratic entropy of a community described by a probability vector and a phylogenetic / functional tree.
}
\usage{
Rao(Ps, Tree, CheckArguments = TRUE)
}
\arguments{
  \item{Ps}{
    A probability vector, summing to 1.
  }
  \item{Tree}{
    An object of class \code{\link{hclust}} or \code{\link{phylog}}. The tree must be ultrametric.
  }
  \item{CheckArguments}{
    Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  Bias correction is not possible when diversity is calculated from probabilities.
  
  Use \code{\link{bcPhyloDiversity}}\code{(Ps, 2, Tree)} to do it.
}
\value{
  A number equal to the calculated diversity.
}
\references{
  Rao, C. R. (1982). Diversity and dissimilarity coefficients: a unified approach. \emph{Theoretical Population Biology} 21(24-43).
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{bcPhyloDiversity}}
}
\examples{
  # Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest
  #      and their taxonomy)
  data(Paracou618)
  # Ps is the vector of probabilities
  Ps <- Paracou618.MC$Ps
  # Calculate the phylogenetic Shannon diversity of the plot
  Rao(Ps, Paracou618.Taxonomy)
}
