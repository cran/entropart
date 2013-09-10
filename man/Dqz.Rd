\name{Dqz}
\alias{Dqz}
\title{
  Similarity-based diversity of a community
}
\description{
  Calculates the diversity of order \eqn{q} of a probability vector according to a similarity matrix.
}
\usage{
Dqz(Ps, q, Z, CheckArguments = TRUE)
}
\arguments{
  \item{Ps}{
    A probability vector, summing to 1.
  }
  \item{q}{
    A number.
  }
  \item{Z}{
    A relatedness matrix, \emph{i.e.} a square matrix whose terms are all positive, strictly positive on the diagonal.
  }
  \item{CheckArguments}{
    Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  Diversity is calculated following Leinster and Cobbold (2012): it is the reciprocal of the (generalized) average (of order \code{q}) of the community species ordinariness.
  
  A similarity matrix is used (as for \code{\link{Dqz}}), not a distance matrix as in Ricotta and Szeidl (2006). See the example.
}
\value{
  A number equal to the calculated diversity.
}
\references{
  Leinster, T. and Cobbold, C. (2012). Measuring diversity: the importance of species similarity. \emph{Ecology} 93(3): 477-489.
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
  DistanceMatrix <- as.matrix(EightSpTree$Wdist^2/2)
  # Similarity can be 1 minus normalized distances between species
  Z <- 1 - DistanceMatrix/max(DistanceMatrix)
  # Calculate diversity of order 2
  Ps <- EightSpAbundance/sum(EightSpAbundance)
  Dqz(Ps, 2, Z)
  # Equal to phylodiversity when q=2
  summary(PhyloDiversity(Ps, 2, EightSpTree))
  # But different for all other q, e.g. 1
  Dqz(Ps, 1, Z)
  summary(PhyloDiversity(Ps, 1, EightSpTree))
}
