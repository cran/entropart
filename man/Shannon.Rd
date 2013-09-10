\name{Shannon}
\alias{bcShannon}
\alias{Shannon}
\title{
  Shannon entropy of a community
}
\description{
  Calculates the Shannon entropy of a probability vector.
}
\usage{
Shannon(Ps, CheckArguments = TRUE)
bcShannon(Ns, Correction = "Best", CheckArguments = TRUE)
}
\arguments{
  \item{Ps}{
    A probability vector, summing to 1.
  }
  \item{Ns}{
    A numeric vector containing species abundances.
  }
  \item{Correction}{
    A string containing one of the possible corrections: \code{"None"} (no correction), \code{"ChaoShen"}, \code{"Grassberger"}, \code{"Grassberger2003"}, \code{"Schurmann"}, \code{"Holste"}, \code{"Bonachela"} or \code{"Best"}, the default value.
  }
  \item{CheckArguments}{
    Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  Bias correction requires the number of individuals to estimate sample \code{\link{Coverage}}. Use \code{bcShannon} and choose the \code{Correction}.
  Correction techniques are from Chao and Shen (2003), Grassberger (1988), Grassberger (2003), Schurmann (2003), Holste et al. (1998), Bonachela et al. (2008).  
  Currently, the \code{"Best"} correction is the max value of \code{"ChaoShen"} and \code{"Grassberger"} (Marcon et al., submitted).
}
\value{
  A number equal to the calculated entropy.
}
\references{
  Bonachela, J. A., Hinrichsen, H. and Munoz, M. A. (2008). Entropy estimates of small data sets. \emph{Journal of Physics A: Mathematical and Theoretical} 41(202001): 1-9.
  
  Chao, A. and Shen, T. J. (2003). Nonparametric estimation of Shannon's index of diversity when there are unseen species in sample. \emph{Environmental and Ecological Statistics} 10(4): 429-443.
  
  Grassberger, P. (1988). Finite sample corrections to entropy and dimension estimates. \emph{Physics Letters A} 128(6-7): 369-373.
  
  Grassberger, P. (2003). Entropy Estimates from Insufficient Samplings. \emph{ArXiv Physics e-prints} 0307138.
  
  Holste, D., Grosse, I. and Herzel, H. (1998). Bayes' estimators of generalized entropies. \emph{Journal of Physics A: Mathematical and General} 31(11): 2551-2566.
  
  Shannon, C. E. (1948). A Mathematical Theory of Communication. \emph{The Bell System Technical Journal} 27: 379-423, 623-656.
  
  Schurmann, T. (2004). Bias analysis in entropy estimation. \emph{Journal of Physics A: Mathematical and Theoretical} 37(27): L295-L301.
  
  Tsallis, C. (1988). Possible generalization of Boltzmann-Gibbs statistics. \emph{Journal of Statistical Physics} 52(1): 479-487.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
    \code{\link{bcShannon}}, \code{\link{Tsallis}}
}
\examples{
  # Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
  data(Paracou618)
  # Ps is the vector of probabilities
  Ps <- Paracou618.MC$Ps
  # Calculate Shannon entropy
  Shannon(Ps)

  # Ns is the vector of abundances of the metacommunity
  Ns <- Paracou618.MC$Ns
  # Calculate an unbiased estimator of Shannon entropy
  bcShannon(Ns)  
}