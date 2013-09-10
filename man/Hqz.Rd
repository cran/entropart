\name{Hqz}
\alias{Hqz}
\title{
  Similarity-based entropy of a community
}
\description{
  Calculates the entropy of order \eqn{q} of a probability vector according to a similarity matrix.
}
\usage{
Hqz(Ps, q, Z, CheckArguments = TRUE)
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
  Entropy is calculated following Leinster and Cobbold (2012) after Ricotta and Szeidl (2006): it is the entropy of order \code{q} of the community, using species ordinariness as the information function.
  
  A similarity matrix is used (as for \code{\link{Dqz}}), not a distance matrix as in Ricotta and Szeidl (2006). See the example.
}
\value{
  A number equal to the calculated entropy.
}
\references{
  Leinster, T. and Cobbold, C. (2012). Measuring diversity: the importance of species similarity. \emph{Ecology} 93(3): 477-489.
  
  Ricotta, C. and Szeidl, L. (2006). Towards a unifying approach to diversity measures: Bridging the gap between the Shannon entropy and Rao's quadratic index. \emph{Theoretical Population Biology} 70(3): 237-243.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{Dqz}}, \code{\link{PhyloEntropy}}
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
  Hqz(Ps, 2, Z)
  # Equal to normalized Rao quadratic entropy when q=2
  Rao(Ps, EightSpTree)/max(DistanceMatrix)
  # But different from PhyloEntropy for all other q, e.g. 1
  Hqz(Ps, 1, Z)
  summary(PhyloEntropy(Ps, 1, EightSpTree))
}
