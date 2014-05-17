\name{Dqz}
\alias{Dqz}
\alias{bcDqz}
\title{
  Similarity-based diversity of a community
}
\description{
  Calculates the diversity of order \eqn{q} of a probability vector according to a similarity matrix.
}
\usage{
Dqz(Ps, q = 1, Z = diag(length(Ps)), CheckArguments = TRUE)
bcDqz(Ns, q = 1, Z = diag(length(Ns)), Correction = "Best", CheckArguments = TRUE)
}
\arguments{
  \item{Ps}{
    A probability vector, summing to 1.
  }
  \item{Ns}{
    A numeric vector containing species abundances.
  }
  \item{q}{
    A number: the order of diversity. Default is 1.
  }
  \item{Z}{
    A relatedness matrix, \emph{i.e.} a square matrix whose terms are all positive, strictly positive on the diagonal. Generally, the matrix is a similarity matrix, \emph{i.e.} the diagonal terms equal 1 and other terms are between 0 and 1. Default is the identity matrix to calculate neutral diversity.
  }
  \item{Correction}{
    A string containing one of the possible corrections: \code{"None"} (no correction), \code{"HorvitzThomson"}, \code{"MarconZhang"} or \code{"Best"}, the default value. The \code{"MarconZhang"} correction assumes a similarity matrix.
  }
  \item{CheckArguments}{
    Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  Diversity is calculated following Leinster and Cobbold (2012): it is the reciprocal of the (generalized) average (of order \code{q}) of the community species ordinariness.
  
  A similarity matrix is used (as for \code{\link{Dqz}}), not a distance matrix as in Ricotta and Szeidl (2006). See the example.

  Bias correction requires the number of individuals. Use \code{bcHqz} and choose the \code{Correction}.
  Correction techniques are from Marcon and Zhang (2014).
  
  Currently, the \code{"Best"} correction is the max value of \code{"HorvitzThomson"} and \code{"MarconZhang"}.
}
\value{
  A number equal to the calculated diversity.
}
\references{
  Leinster, T. and Cobbold, C. (2012). Measuring diversity: the importance of species similarity. \emph{Ecology} 93(3): 477-489.

  Marcon, E. and Zhang, Z. (2014). The decomposition of similarity-based diversity and its bias correction. \emph{HAL} hal-00989454(version 1).
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{Hqz}}, \code{\link{PhyloDiversity}}
}
\examples{
  # Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
  data(Paracou618)
  # Prepare the similarity matrix
  DistanceMatrix <- as.matrix(Paracou618.dist)
  # Similarity can be 1 minus normalized distances between species
  Z <- 1 - DistanceMatrix/max(DistanceMatrix)
  # Calculate diversity of order 2
  bcDqz(Paracou618.MC$Ns, 2, Z)
}
