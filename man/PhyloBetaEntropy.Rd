\name{PhyloBetaEntropy}
\alias{PhyloBetaEntropy}
\alias{bcPhyloBetaEntropy}
\title{
  Phylogenetic Beta Entropy of a community
}
\description{
  Calculates the phylogenetic beta entropy of order \eqn{q} of a a community belonging to a metacommunity.
}
\usage{
PhyloBetaEntropy(Ps, Pexp, q, Tree, Normalize = TRUE, CheckArguments = TRUE)
bcPhyloBetaEntropy(Ns, Nexp, q, Tree, Normalize = TRUE, Correction = "Best", 
  CheckArguments = TRUE)
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
  \item{Tree}{
    An object of class \code{\link{hclust}} or \code{\link{phylog}}. The tree must be ultrametric.
  }
  \item{Normalize}{
    If \code{TRUE} (default), the entropy returned by the function is normalized by the height of the tree (it is the weighted average value of the entropy in each slice).\cr
    If \code{FALSE}, it is the unnormalized weighted sum of the results.
  }
  \item{Correction}{
    A string containing one of the possible corrections: currently, only \code{"ChaoShen"}. \code{"Best"} is the default value, it is equivalent to \code{"ChaoShen"}.
  }
  \item{CheckArguments}{
    Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  The phylogenetic entropy is the generalization of HCDT entropy to unequal species distances (Pavoine et al., 2009).
  
  Calculation relies on \code{\link{TsallisBeta}} and \code{\link{PhyloApply}}.
  
  Bias correction requires the number of individuals to estimate sample \code{\link{Coverage}}. Use \code{bcPhyloBetaEntropy} and choose the \code{Correction}.
  
  Note that beta entropy value is related to alpha entropy (if \eqn{q} is not 1) and cannot be compared accross communities (Jost, 2007). Beta entropy of a community is not meaningful in general, do rather calculate the \code{\link{PhyloDiversity}} of the metacommunity.
}
\value{
  A \code{\link{PhyloEntropy}} object containing entropy values at each cut of the tree.
}
\references{
    Jost (2007), Partitioning diversity into independent alpha and beta components. \emph{Ecology} 88(10): 2427-2439.

  Marcon, E., Herault, B. (in prep). Partitioning phylodiversity. \emph{Methods in Ecology and Evolution}.
  
  Pavoine, S., Love, M. S. and Bonsall, M. B. (2009). Hierarchical partitioning of evolutionary and ecological patterns in the organization of phylogenetically-structured species assemblages: Application to rockfish (genus: Sebastes) in the Southern California Bight. \emph{Ecology Letters} 12(9): 898-908.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{TsallisBeta}}, \code{\link{bcPhyloBetaEntropy}}, \code{\link{PhyloDiversity}}
}
\examples{
  # Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest
  #      and their taxonomy)
  data(Paracou618)
  # Ps is the vector of probabilities
  Ps <- Paracou618.MC$Ps
  # Probability distribution of the first plot
  Ps1 <- Paracou618.MC$Psi[, 1]
  # Calculate the phylogenetic Shannon beta entropy of the plot
  summary(PhyloBetaEntropy(Ps1, Ps, 1, Paracou618.Taxonomy) -> e)
  plot(e)
  
  # Ns is the vector of abundances of the metacommunity
  Ns <- Paracou618.MC$Ns
  # Abundances in the first plot
  Ns1 <- Paracou618.MC$Nsi[, 1]
  # Calculate the phylogenetic Shannon beta entropy of the plot
  summary(bcPhyloBetaEntropy(Ns1, Ns, 1, Paracou618.Taxonomy, Correction = "Best") -> e)
  plot(e)
}
