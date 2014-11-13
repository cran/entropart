\name{ShannonBeta}
\alias{ShannonBeta}
\alias{bcShannonBeta}
\title{
  Shannon beta entropy of a community
}
\description{
  Calculates the Shannon beta entropy of a community belonging to a metacommunity.
}
\usage{
ShannonBeta(Ps, Pexp, CheckArguments = TRUE)
bcShannonBeta(Ns, Nexp, Correction = "Best", CheckArguments = TRUE)
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
   A string containing one of the possible corrections: currently, \code{"ChaoShen"} (Marcon \emph{et al.}, 2012) equivalent to \code{"Best"}, and \code{"ZhangGrabchak"} (Zhang and Grabchak, 2014).
  }
  \item{CheckArguments}{
    Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  The derivation of Shannon beta entropy can be found in Marcon \emph{et al.} (2012).
  
  Bias correction requires the number of individuals to estimate sample \code{\link{Coverage}}. Use \code{bcShannonBeta} and choose the \code{Correction}.

}
\value{
  A number equal to the calculated entropy.
}
\references{
  Marcon, E., Herault, B., Baraloto, C. and Lang, G. (2012). The Decomposition of Shannon's Entropy and a Confidence Interval for Beta Diversity. \emph{Oikos} 121(4): 516-522.
  
  Zhang, Z. and Grabchak M. (2014). Nonparametric Estimation of Kullback-Leibler Divergence. \emph{Neural computation} in press.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{bcShannonBeta}}
}
\examples{
  # Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
  data(Paracou618)
  # Ps is the vector of probabilities
  Ps <- Paracou618.MC$Ps
  # Probability distribution of the first plot
  Ps1 <- Paracou618.MC$Psi[, 1]
  # Shannon beta entropy of the plot
  ShannonBeta(Ps1, Ps)

  # Ns is the vector of abundances of the metacommunity
  Ns <- Paracou618.MC$Ns
  # Abundances in the first plot
  Ns1 <- Paracou618.MC$Nsi[, 1]
  # Reduced-bias estimator of Shannon beta entropy of the plot
  bcShannonBeta(Ns1, Ns)
}
