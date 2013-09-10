\name{Diversity}
\alias{Diversity}
\alias{bcDiversity}
\title{
  Generalized diversity of a community
}
\description{
  Calculates the generalized diversity of order \eqn{q} of a probability vector.
}
\usage{
Diversity(Ps, q, CheckArguments = TRUE)
bcDiversity(Ns, q, Correction = "Best", CheckArguments = TRUE)
}
\arguments{
  \item{Ps}{
    A probability vector, summing to 1.
  }
  \item{Ns}{
    A numeric vector containing species abundances.
  }
  \item{q}{
    A number.
  }
  \item{Correction}{
    A string containing one of the possible corrections: \code{"None"} (no correction), \code{"ChaoShen"}, \code{"Grassberger"}, \code{"Holste"}, \code{"Bonachela"} or \code{"Best"}, the default value.
  }
  \item{CheckArguments}{
    Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  \code{Diversity} calls \code{\link{Tsallis}} to calculate entropy and transforms it into diversity by calculating its deformed exponential.
  
  Bias correction requires the number of individuals to estimate sample \code{\link{Coverage}}. Use \code{bcDiversity} and choose the \code{Correction}.
}
\value{
  A number equal to the calculated diversity.
}
\references{
  Marcon, E., Scotti, I., Herault, B., Rossi, V. and Lang, G. (in revision). Generalization of the partitioning of Shannon diversity. \emph{PLOS One}.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{Tsallis}}, \code{\link{expq}}
}
\examples{
  # Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
  data(Paracou618)
  # Ns is the total number of trees per species
  Ns <- Paracou618.MC$Ns
  # Whittaker plot
  plot(sort(Ns[Ns != 0], decr = TRUE), log = "y", type = "b",
    main = "Whittaker plot", xlab = "Rank", ylab = "Ns")
  # Calculate diversity of order 1, i.e. Shannon's diversity
  Diversity(Paracou618.MC$Ps, 1)
  # Calculate it with estimation bias correction
  bcDiversity(Paracou618.MC$Ns, 1)
}
