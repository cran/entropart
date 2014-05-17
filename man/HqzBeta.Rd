\name{HqzBeta}
\alias{HqzBeta}
\alias{bcHqzBeta}
\title{
  Similarity-based beta entropy of a community
}
\description{
  Calculates the similarity-based beta entropy of order \eqn{q} of a community belonging to a metacommunity.
}
\usage{
HqzBeta(Ps, Pexp = NULL, q = 1, Z = diag(length(Ps)), CheckArguments = TRUE) 
bcHqzBeta(Ns, Nexp = NULL, q = 1, Z = diag(length(Ns)), Correction = "Best",
          CheckArguments = TRUE)
}
\arguments{
  \item{Ps}{
  The probability vector of species of the community.
  }
  \item{Pexp}{
  The probability vector of species of the metacommunity.
  }
  \item{Ns}{
    A numeric vector containing species abundances of the community.
  }
  \item{Nexp}{
    A numeric vector containing species abundances of the metacommunity.
  }
  \item{q}{
  A number, the order of diversity. Default is 1.
  }
  \item{Z}{
    A relatedness matrix, \emph{i.e.} a square matrix whose terms are all positive, strictly positive on the diagonal. Generally, the matrix is a similarity matrix, \emph{i.e.} the diagonal terms equal 1 and other terms are between 0 and 1. Default is the identity matrix to calculate neutral entropy.
  }
  \item{Correction}{
    A string containing one of the possible corrections: currently, no correction is available so \code{"Best"},  the default value, is equivalent to \code{"None"}.
  }
  \item{CheckArguments}{
    Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  The derivation of similarity-based beta entropy can be found in Marcon and Zhang (2014).
  
  Bias correction requires the number of individuals.

  Note that beta entropy value is related to alpha entropy (if \eqn{q} is not 1) and cannot be compared accross communities (Jost, 2007). Beta entropy of a community is not meaningful in general, do rather calculate the \code{\link{BetaDiversity}} of the metacommunity.
}
\value{
  A number equal to the calculated entropy.
}
\references{
  Jost (2007), Partitioning diversity into independent alpha and beta components. \emph{Ecology} 88(10): 2427-2439.

  Marcon, E. and Zhang, Z. (2014). The decomposition of similarity-based diversity and its bias correction. \emph{HAL} hal-00989454(version 1).
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\examples{
  # Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
  data(Paracou618)
  # Ps is the vector of probabilities
  Ps <- Paracou618.MC$Ps
  # Probability distribution of the first plot
  Ps1 <- Paracou618.MC$Psi[, 1]
  # Prepare the similarity matrix
  DistanceMatrix <- as.matrix(Paracou618.dist)
  # Similarity can be 1 minus normalized distances between species
  Z <- 1 - DistanceMatrix/max(DistanceMatrix)
  # Divergence of order 2 between plot 1 and the whole forest
  HqzBeta(Ps1, Ps, 2, Z)
}
