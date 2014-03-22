\name{Richness}
\alias{Richness}
\title{
  Simpson entropy of a community
}
\description{
  Calculates the richness of a community (described by a vector of abundances or probabilities).
}
\usage{
Richness(Ns, CheckArguments = TRUE)
}
\arguments{
  \item{Ns}{
    A numeric vector. Contains either abundances or probabilities.
  }
  \item{CheckArguments}{
    Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  Richness simply returns the number of strictly positive elements of the vector.
}
\value{
  A number equal to the calculated entropy.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{Tsallis}}
}
\examples{
  # Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
  data(Paracou618)
  # Richness of the first community of Paracou618.MC
  Richness(Paracou618.MC$Psi[, 1]) 
}