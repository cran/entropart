\name{Preprocess.Tree}
\alias{Preprocess.Tree}
\title{
  Preprocessing of a phylogenetic tree
}
\description{
  Calculates cuts and intervals of a phylogenetic tree and make it available both in \code{\link{hclust}} and \code{\link{phylog}} formats.
}
\usage{
Preprocess.Tree(Tree)
}
\arguments{
  \item{Tree}{
    An object of class \code{\link{hclust}} or \code{\link{phylog}}. The tree must be ultrametric.
  }
}
\details{
  This is an internal function. It is used by \code{\link{PhyloApply}} and to obtain exact values of intervals (functions such as \code{\link{PhyloDiversity}} return them as names, with low precision).
}
\value{
  An object of class \code{\link{PPtree}}.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\keyword{internal}
