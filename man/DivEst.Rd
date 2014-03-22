\name{DivEst}
\alias{DivEst}
\alias{plot.DivEst}
\alias{summary.DivEst}
\title{
  Diversity Estimation of a metacommunity
}
\description{
  Estimates diversity of a metacommunity.
}
\usage{
DivEst(q = 0, MC, Biased = TRUE, Correction = "Best", Tree = NULL, 
  Normalize = TRUE, Simulations = 100, CheckArguments = TRUE)
\method{plot}{DivEst}(x, \dots, main = NULL, Which = "All")
\method{summary}{DivEst}(object, \dots)
}
\arguments{
  \item{q}{
    A number: the order of diversity.
  }
  \item{MC}{
    A \code{\link{MetaCommunity}} object.
  }
  \item{Biased}{
    Logical; if \code{FALSE}, a bias correction is appplied.
  }
  \item{Correction}{
    A string containing one of the possible corrections. The correction must be accepted by \code{\link{DivPart}}. \code{"Best"} is the default value.
  }
  \item{Tree}{
    An object of class \code{\link{hclust}} or \code{\link{phylog}}. The tree must be ultrametric.
  }
  \item{Normalize}{
    If \code{TRUE} (default), diversity is not affected by the height of the tree..\cr
    If \code{FALSE}, diversity is proportional to the height of the tree.
  }
  \item{Simulations}{
    The number of simulations to build confidence intervals.
  }
  \item{CheckArguments}{
    Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
  \item{x}{
    An object to be tested or plotted.
  }
  \item{main}{
    The title of the plot.
  }
  \item{Which}{
    May be \code{"Alpha"}, \code{"Beta"} or \code{"Gamma"} to respectively plot the metacommunity's alpha, beta or gamma diversity. If \code{"All"} (default), all three plots are shown.
    }
  \item{object}{
    A \code{MCdiversity} object to be summarized.
  }
  \item{\dots}{
    Additional arguments to be passed to the generic methods.
  }
}
\details{
  \code{Divest} estimates the diversity of the metacommunity and partitions it into alpha and beta components.
  
  If \code{Tree} is provided, the phylogenetic diversity is calculated.
  
  Confidence intervals are calculated by Monte-Carlo simulations, drawing simulated communities from a multinomial law following observed frequencies (Marcon et al, 2012 ; submitted)
}
\value{
  A \code{Divest} object which is a \code{\link{DivPart}} object with an additional item in its list:
    \item{SimulatedDiversity}{A matrix containing the simulated values of alpha, beta and gamma diversity.}

  \code{Divest} objects can be summarized and plotted.
}
\references{
  Marcon, E., Herault, B., Baraloto, C. and Lang, G. (2012). The Decomposition of Shannon's Entropy and a Confidence Interval for Beta Diversity. \emph{Oikos} 121(4): 516-522.
  
  Marcon, E., Scotti, I., Herault, B., Rossi, V. and Lang, G. (2014). Generalization of the partitioning of Shannon diversity. \emph{PLOS One} 9(3): e90289.

  Marcon, E., Herault, B. (2014). Decomposing  Phylodiversity. \emph{HAL} hal-00946177(version 1).
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf> and Bruno Herault <Bruno.Herault@ecofog.gf>
}
\seealso{
  \code{\link{DivPart}}
}
\examples{
  # Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
  data(Paracou618)
  # Estimate Shannon diversity.
  Estimation <- DivEst(q = 1, Paracou618.MC, Biased = FALSE, Correction = "Best", 
    Simulations = 1000)
  plot(Estimation)
  summary(Estimation)
}
