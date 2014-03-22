\name{Simpson}
\alias{Simpson}
\alias{bcSimpson}
\title{
  Simpson entropy of a community
}
\description{
  Calculates the Simpson entropy of a probability vector.
}
\usage{
Simpson(Ps, CheckArguments = TRUE)
bcSimpson(Ns, Correction = "Lande", CheckArguments = TRUE)
}
\arguments{
  \item{Ps}{
    A probability vector, summing to 1.
  }
  \item{Ns}{
    A numeric vector containing species abundances.
  }
  \item{Correction}{
    A string containing one of the possible corrections accepted by \code{\link{bcTsallis}} or \code{"Lande"}, the default value.
  }
  \item{CheckArguments}{
    Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  Lande's correction has been derived (Lande, 1996; Good, 1953) especially for Simpson entropy, while other corrections are for generalized Tsallis entropy. It is identical to the unbiased estimator proposed by Simpson, although arguments were different.
  Bias correction requires the number of individuals to estimate sample \code{\link{Coverage}}. Use \code{bcSimpson} and choose the \code{Correction}.
}
\value{
  A number equal to the calculated entropy.
}
\references{
  Good, I. J. (1953). On the Population Frequency of Species and the Estimation of Population Parameters. \emph{Biometrika} 40(3/4): 237-264.

  Lande, R. (1996). Statistics and partitioning of species diversity, and similarity among multiple communities. \emph{Oikos} 76: 5-13.

  Simpson, E. H. (1949). Measurement of diversity. \emph{Nature} 163(4148): 688.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{Tsallis}}, \code{\link{bcSimpson}}
}
\examples{
  # Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
  data(Paracou618)
  # Ps is the vector of probabilities
  Ps <- Paracou618.MC$Ps
  # Calculate Simpson entropy
  Simpson(Ps)

  # Ns is the vector of abundances of the metacommunity
  Ns <- Paracou618.MC$Ns
  # Calculate an unbiased estimator of Simpson's index of diversity
  bcSimpson(Ns)  
}