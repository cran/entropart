\name{DivProfile}
\alias{DivProfile}
\alias{plot.DivProfile}
\alias{summary.DivProfile}
\title{
  Diversity Profile of a metacommunity
}
\description{
  Calculate the diversity profiles (alpha, beta, gamma) of a metacommunity.
}
\usage{
DivProfile(q.seq = seq(0, 2, 0.1), MC, Biased = TRUE, Correction = "Best", 
  Tree = NULL, Normalize = TRUE, CheckArguments = TRUE)
\method{plot}{DivProfile}(x, \dots, main = NULL, xlab = "Order of Diversity",
  ylab = NULL, Which = "All")
\method{summary}{DivProfile}(object, \dots)
}
\arguments{
  \item{q.seq}{
    A numeric vector.
  }
  \item{MC}{
    A \code{\link{MetaCommunity}} object.
  }
  \item{Biased}{
    Logical; if \code{FALSE}, a bias correction is appplied.
  }
  \item{Correction}{
    A string containing one of the possible corrections.
	
    The correction must be accepted by \code{\link{AlphaEntropy}}, \code{\link{BetaEntropy}} and \code{\link{GammaEntropy}}. \code{"Best"} is the default value.
  }
  \item{Tree}{
    An object of class \code{\link{hclust}} or \code{\link{phylog}}. The tree must be ultrametric.
  }
  \item{Normalize}{
    If \code{TRUE} (default), diversity is not affected by the height of the tree.\cr
    If \code{FALSE}, diversity is proportional to the height of the tree.
  }
  \item{CheckArguments}{
    Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
  \item{x}{
    An object to be tested or plotted.
  }
  \item{main}{
    The main title of the plot. Ignored if \code{Which = "All"}.
  }
  \item{xlab}{
    The x axis label of the plots.
  }
  \item{ylab}{
    The y axis label of the plot. Ignored if \code{Which = "All"}.
  }
  \item{Which}{
    May be \code{"Communities"}, \code{"Alpha"}, \code{"Beta"} or \code{"Gamma"} to respectively plot the alpha diversity of communities or the metacommunity's alpha, beta or gamma diversity. If \code{"All"} (default), all four plots are shown.
  }
  \item{object}{
    A \code{MCdiversity} object to be summarized.
  }
  \item{\dots}{
    Additional arguments to be passed to the generic methods.
  }
}
\details{
  If \code{Tree} is provided, the phylogenetic diversity is calculated.
  
  \code{DivPart} partitions the diversity of the metacommunity into alpha and beta components. It supports estimation-bias correction (Marcon et al., submitted).
  
  If \code{Tree} is provided, the phylogenetic diversity is calculated.
  
  Beta diversity/entropy is calculated from Gamma and Alpha when bias correction is required, so community values are not available.
}
\value{
  A \code{DivProfile} object. It is a list:
    \item{MetaCommunity}{The name of the \code{\link{MetaCommunity}} object containing inventory data.}
    \item{Order}{A vector containing the values of \code{q}.}
    \item{Biased}{Logical. If \code{FALSE}, bias corrected values of diversity have been computed.}
    \item{Correction}{The estimation bias correction used to calculate diversity.}
    \item{Tree}{The phylogenetic or functional tree used to calculate phylodiversity.}
    \item{Normalized}{Logical. Indicates whether phylodiversity is normalized or proportional to the height of the tree.}
    \item{CommunityAlphaDiversities}{A matrix containing the alpha diversity of each community.}
    \item{TotalAlphaDiversity}{A vector containing the alpha diversity of communities for each order.}
    \item{BetaDiversity}{A vector containing the beta diversity of communities for each order.}
    \item{GammaDiversity}{A vector containing the gamma diversity of the metacommunity for each order.}
    \item{CommunityAlphaEntropies}{A matrix containing the alpha entropy of each community.}
    \item{TotalAlphaEntropy}{A vector containing the alpha entropy of communities for each order.}
    \item{BetaEntropy}{A vector containing the beta entropy of communities for each order.}
    \item{GammaEntropy}{A vector containing the gamma entropy of the metacommunity for each order.}
  
  \code{DivProfile} objects can be summarized and plotted.
}
\references{
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
  # Estimate diversity.
  Profile <- DivProfile(q.seq = seq(0, 2, 0.1), Paracou618.MC, Biased = FALSE)
  plot(Profile)
  summary(Profile)
}
