\name{CommunityProfile}
\alias{CommunityProfile}
\title{
  Diversity or Entropy Profile of a community
}
\description{
  Calculates the diversity or entropy profile of a community, applying a community function to a vector of orders.
}
\usage{
CommunityProfile(FUN, NorP, q.seq, ..., CheckArguments = TRUE)
}
\arguments{
  \item{FUN}{
    The function to be applied to each value of \eqn{q}. 
  }
  \item{NorP}{
    A numeric vector or a two-column matrix. Contains either abundances or probabilities. Two-column matrices should contain the observed abundances (or probabilities) in the first column and the expected ones in the second column, to allow using beta diversity functions.
  }
  \item{q.seq}{
    A numeric vector.
  }
  \item{...}{
    Additional arguments to be passed to \code{FUN}.
  }
  \item{CheckArguments}{
    Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  This function is used to calculate diversity or entropy profiles based on community functions such as \code{\link{Tsallis}} or \code{\link{ChaoPD}}. The first two arguments of the function must be a probability or abundance vector or a two-column matrix (for beta diversity functions) and a number (\eqn{q}).
  Additional arguments cannot be checked. Unexpected results may be returned if \code{FUN} is not used properly.
}
\value{
  A list:
    \item{x}{The order \eqn{q} values}
    \item{y}{The entropy or diversity values returned by \code{FUN}}
  The list can be plotted.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\examples{
  # Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
  data(Paracou618)
  # Estimate diversity.
  Profile <- CommunityProfile(bcDiversity, Paracou618.MC$Ns, seq(0, 2, 0.1))
  plot(Profile, type="l", main="Paracou Plots Diversity", xlab="q", ylab="Diversity")
  # Compare to non bias-corrected diversity
  lines(CommunityProfile(Diversity, Paracou618.MC$Ps, seq(0, 2, 0.1)), lty=3)
  legend("topright", c("Biased", "Bias Corrected"), lty=c(1,3), inset=0.01)
}
