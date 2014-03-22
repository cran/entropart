\name{plot.PhyloValue}
\alias{plot.PhyloValue}
\title{
  Plot phylogenetic entropy of diversity.
}
\description{
  Plots entropy or diversity against the height of the phylogenetic or functional tree.
}
\usage{
\method{plot}{PhyloValue}(x, xlab = expression(italic("T")), ylab = NULL, main = NULL, \dots)
}
\arguments{
  \item{x}{
    An object of class \code{PhyloValue}, including \code{\link{PhyloDiversity}} and \code{\link{PhyloEntropy}} objects.
  }
  \item{xlab}{
    The X axis label, "T" by default for Time.
  }
  \item{ylab}{
    The Y axis label. if \code{NULL} (by default), "Entropy" or "Diversity" or nothing is chosen according to the object class.
  }
  \item{main}{
    The main title of the plot. if \code{NULL} (by default), a default value is used.
  }
  \item{\dots}{
    Additional arguments to be passed to \code{\link{plot}}.
  }
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
