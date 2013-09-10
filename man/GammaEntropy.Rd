\name{GammaEntropy}
\alias{GammaEntropy}
\title{
  Unbiased gamma entropy of a metacommunity
}
\description{
  Calculates the unbiased Tsallis entropy of order \eqn{q} of a metacommunity.
}
\usage{
GammaEntropy(MC, q, Correction = "Best", Tree = NULL, Normalize = TRUE, 
  CheckArguments = TRUE)
}
\arguments{
  \item{MC}{
    A \code{\link{MetaCommunity}} object.
  }
  \item{q}{
    A number: the order of diversity.
  }
  \item{Correction}{
    A string containing one of the possible corrections: \code{"ChaoShen"}, \code{"Grassberger"}, \code{"Holste"}, \code{"Bonachela"} or \code{"Best"}, the default value.
  }
  \item{Tree}{
    An object of class \code{\link{hclust}} or \code{\link{phylog}}. The tree must be ultrametric.
  }
  \item{Normalize}{
    If \code{TRUE} (default), the entropy returned by the function is normalized by the height of the tree (it is the weighted average value of the entropy in each slice).\cr
    If \code{FALSE}, it is the unnormalized weighted sum of the results.
  }
  \item{CheckArguments}{
    Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  If \code{Tree} is \code{NULL}, then generalized entropy is calculated by \code{\link{bcTsallis}}, else phylogenetic entropy is calculated by \code{\link{bcPhyloEntropy}}.
  
  Correction techniques are from Chao and Shen (2003), Grassberger (1988), Holste et al. (1998), Bonachela et al. (2008).
  Currently, the \code{"Best"} correction is the max value of \code{"ChaoShen"} and \code{"Grassberger"} (Marcon et al., submitted).
}
\value{
  A number equal to the calculated entropy.
}
\references{
  Bonachela, J. A., Hinrichsen, H. and Munoz, M. A. (2008). Entropy estimates of small data sets. \emph{Journal of Physics A: Mathematical and Theoretical} 41(202001): 1-9.
  
  Chao, A. and Shen, T. J. (2003). Nonparametric estimation of Shannon's index of diversity when there are unseen species in sample. \emph{Environmental and Ecological Statistics} 10(4): 429-443.
  
  Grassberger, P. (1988). Finite sample corrections to entropy and dimension estimates. \emph{Physics Letters A} 128(6-7): 369-373.
  
  Holste, D., Grosse, I. and Herzel, H. (1998). Bayes' estimators of generalized entropies. \emph{Journal of Physics A: Mathematical and General} 31(11): 2551-2566.
  
  Marcon, E., Scotti, I., Herault, B., Rossi, V. and Lang, G. (in revision). Generalization of the partitioning of Shannon diversity. \emph{PLOS One}.
  
  Marcon, E., Herault, B. (in prep). Partitioning phylodiversity. \emph{Methods in Ecology and Evolution}.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{bcTsallis}}, \code{\link{bcPhyloEntropy}}
}
\examples{
  # Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
  data(Paracou618)
  # Calculate Simpson gamma entropy
  GammaEntropy(Paracou618.MC, 2)
  # Compare without correction
  GammaEntropy(Paracou618.MC, 2, Correction = "None")
  # Estimate phylogenetic Simpson gamma entropy
  GammaEntropy(Paracou618.MC, 2, Tree = Paracou618.Taxonomy)   
}
