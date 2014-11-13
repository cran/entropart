\name{Rao}
\alias{Rao}
\alias{bcRao}
\title{
  Rao Quadratic Entropy of a Community
}
\description{
  Calculates Rao's quadratic entropy of a community described by a probability vector and a phylogenetic / functional tree.
}
\usage{
Rao(Ps, Tree, CheckArguments = TRUE)
bcRao(Ns, Tree, Correction="Lande", CheckArguments = TRUE)
}
\arguments{
  \item{Ps}{
    A probability vector, summing to 1.
  }
  \item{Ns}{
    A numeric vector containing species abundances.
  }
  \item{Tree}{
    An object of class \code{\link{hclust}} or \code{\link{phylog}}. The tree must be ultrametric.
  }
  \item{Correction}{
    A string containing one of the possible corrections accepted by \code{\link{bcTsallis}} or \code{"Lande"}, the default value (equivalent to \code{"Best"}).
  }
  \item{CheckArguments}{
    Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
    Bias correction requires the number of individuals. Use \code{bcRao} and choose the \code{Correction}.
    
    The unbiased estimator of Rao's entropy is identical to that of Simpson's entropy because Rao's entropy is a linear sum of Simson entropies, all of them calculated from the same number of individuals (Marcon and Herault, 2014). It equals the plug-in etimator multiplied by n/(n-1) where n is the total number of individuals.

}
\value{
  A number equal to the calculated diversity.
}
\references{
  Marcon, E., Herault, B. (2014). Decomposing  Phylodiversity. \emph{HAL} hal-00946177(version 1).

  Rao, C. R. (1982). Diversity and dissimilarity coefficients: a unified approach. \emph{Theoretical Population Biology} 21(24-43).
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{bcPhyloDiversity}}
}
\examples{
  # Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest
  #      and their taxonomy)
  data(Paracou618)
  # Ps is the vector of probabilities
  Ps <- Paracou618.MC$Ps
  # Calculate the phylogenetic Shannon diversity of the plot
  Rao(Ps, Paracou618.Taxonomy)
}
