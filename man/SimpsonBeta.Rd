\name{SimpsonBeta}
\alias{SimpsonBeta}
\alias{bcSimpsonBeta}
\title{
  Simpson beta entropy of a community
}
\description{
  Calculates the Simpson beta entropy of a community belonging to a metacommunity.
}
\usage{
SimpsonBeta(Ps, Pexp, CheckArguments = TRUE)
bcSimpsonBeta(Ns, Nexp, Correction = "Best", CheckArguments = TRUE)
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
  \item{Correction}{
    A string containing one of the possible corrections: currently, only \code{"ChaoShen"}, identical to \code{"Best"}.
  }
  \item{CheckArguments}{
    Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  The derivation of Tsallis beta entropy (Simpson is Tsallis of order 2) can be found in Marcon et al. (in revision).
  
  Bias correction requires the number of individuals to estimate sample \code{\link{Coverage}}. Use \code{bcSimpsonBeta} and choose the \code{Correction}.
  
  Note that Simpson beta entropy value is related to Simpson alpha entropy value and cannot be compared accross communities (Jost, 2007). Beta entropy of a community is not meaningful in general, do rather calculate the \code{\link{BetaDiversity}} of order 2 of the metacommunity.
}
\value{
  A number equal to the calculated entropy.
}
\references{
  Jost (2007), Partitioning diversity into independent alpha and beta components. \emph{Ecology} 88(10): 2427-2439.
  
  Marcon, E., Scotti, I., Herault, B., Rossi, V. and Lang, G. (in revision). Generalization of the partitioning of Shannon diversity. \emph{PLOS One}.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{Simpson}}, \code{\link{bcSimpsonBeta}}, \code{\link{BetaDiversity}}
}
\examples{
  # Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
  data(Paracou618)
  # Ps is the vector of probabilities
  Ps <- Paracou618.MC$Ps
  # Probability distribution of the first plot
  Ps1 <- Paracou618.MC$Psi[, 1]
  # Shannon beta entropy of the plot
  SimpsonBeta(Ps1, Ps)
  # Transform into diversity
  expq(SimpsonBeta(Ps1, Ps)/(1-Simpson(Ps1)), 2)
  
  # Ns is the vector of abundances of the metacommunity
  Ns <- Paracou618.MC$Ns
  # Abundances in the first plot
  Ns1 <- Paracou618.MC$Nsi[, 1]
  # Unbiased Shannon beta entropy of the plot
  bcSimpsonBeta(Ns1, Ns)
}
