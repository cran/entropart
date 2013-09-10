\name{Tsallis}
\alias{Tsallis}
\alias{bcTsallis}
\title{
  Tsallis Entropy of a community
}
\description{
  Calculates the HCDT, also known as Tsallis entropy of order \eqn{q} of a probability vector.
}
\usage{
Tsallis(Ps, q, CheckArguments = TRUE)
bcTsallis(Ns, q, Correction = "Best", CheckArguments = TRUE)
}
\arguments{
  \item{Ps}{
    A probability vector, summing to 1.
  }
  \item{Ns}{
    A numeric vector containing species abundances.
  }
  \item{q}{
    A number: the order of entropy. Some corrections allow only a positive number.
  }
  \item{Correction}{
    A string containing one of the possible corrections: \code{"None"} (no correction), \code{"ChaoShen"}, \code{"Grassberger"}, \code{"Holste"}, \code{"Bonachela"} or \code{"Best"}, the default value.
  }
  \item{CheckArguments}{
    Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  Tsallis (Havrda and Charvat, 1967; Daroczy, 1970; Tsallis, 1988) generalized entropy is a generalized measure of diversity (Jost, 2006).
  
  Bias correction requires the number of individuals to estimate sample \code{\link{Coverage}}. Use \code{bcTsallis} and choose the \code{Correction}.
  Correction techniques are from Chao and Shen (2003), Grassberger (1988), Holste et al. (1998), Bonachela et al. (2008).
  
  Currently, the \code{"Best"} correction is the max value of \code{"ChaoShen"} and \code{"Grassberger"} (Marcon et al., submitted).
}
\value{
  A number equal to the calculated entropy.
}
\references{
  Havrda, J. and Charvat, F. (1967). Quantification method of classification processes. Concept of structural a-entropy. \emph{Kybernetika} 3(1): 30-35.
  
  Daroczy, Z. (1970). Generalized information functions. \emph{Information and Control} 16(1): 36-51.
  
  Tsallis, C. (1988). Possible generalization of Boltzmann-Gibbs statistics. \emph{Journal of Statistical Physics} 52(1): 479-487.
  
  Jost, L. (2006). Entropy and diversity. \emph{Oikos} 113(2): 363-375.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{bcTsallis}}
}
\examples{
  # Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
  data(Paracou618)
  # Ps is the vector of probabilities
  Ps <- Paracou618.MC$Ps
  # Calculate Tsallis entropy of order 2, i.e. Simpson's index of diversity
  Tsallis(Ps, 2) 

  # Ns is the vector of abundances of the metacommunity
  Ns <- Paracou618.MC$Ns
  # Calculate an unbiased estimator of Tsallis entropy of order 2, i.e. Simpson's index
  bcTsallis(Ns, 2)
}
