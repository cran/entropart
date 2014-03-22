\name{TsallisBeta}
\alias{TsallisBeta}
\alias{bcTsallisBeta}
\title{
  Tsallis beta entropy of a community
}
\description{
  Calculates the Tsallis beta entropy of order \eqn{q} of a community belonging to a metacommunity.
}
\usage{
TsallisBeta(Ps, Pexp = NULL, q, CheckArguments = TRUE)
bcTsallisBeta(Ns, Nexp = NULL, q, Correction = "Best", CheckArguments = TRUE)
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
  A number.
  }
  \item{Correction}{
    A string containing one of the possible corrections: currently, only \code{"ChaoShen"}. \code{"Best"} is the default value, it is equivalent to \code{"ChaoShen"}.
  }
  \item{CheckArguments}{
    Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  The derivation of Tsallis beta entropy can be found in Marcon et al. (2014).
  
  Bias correction requires the number of individuals to estimate sample \code{\link{Coverage}}. Use \code{bcTsallisBeta} and choose the \code{Correction}.

  Note that beta entropy value is related to alpha entropy (if \eqn{q} is not 1) and cannot be compared accross communities (Jost, 2007). Beta entropy of a community is not meaningful in general, do rather calculate the \code{\link{BetaDiversity}} of the metacommunity.
}
\value{
  A number equal to the calculated entropy.
}
\references{
    Jost (2007), Partitioning diversity into independent alpha and beta components. \emph{Ecology} 88(10): 2427-2439.

  Marcon, E., Scotti, I., Herault, B., Rossi, V. and Lang, G. (2014). Generalization of the partitioning of Shannon diversity. \emph{PLOS One} 9(3): e90289.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{bcTsallisBeta}}
}
\examples{
  # Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
  data(Paracou618)
  # Ps is the vector of probabilities
  Ps <- Paracou618.MC$Ps
  # Probability distribution of the first plot
  Ps1 <- Paracou618.MC$Psi[, 1]
  # Divergence of order 2 between plot 1 and the whole forest
  TsallisBeta(Ps1, Ps, 2)

  # Ns is the vector of abundances of the metacommunity
  Ns <- Paracou618.MC$Ns
  # Abundances in the first plot
  Ns1 <- Paracou618.MC$Nsi[, 1]
  # Divergence of order 2 between plot 1 and the whole forest, with bias correction
  bcTsallisBeta(Ns1, Ns, 2)
}
