\name{Originality.Species}
\alias{Originality.Species}
\title{
  Originality of species
}
\description{
  Calculates the originality of species according to the order of diversity.
}
\usage{
Originality.Species(Tree, q = 2, CheckArguments = TRUE)
}
\arguments{
  \item{Tree}{
    An object of class \code{\link{hclust}} or \code{\link{phylog}}. The tree must be ultrametric.
  }
  \item{q}{
    A number: the order of diversity. Default is 2 for Rao quadratic entropy.
  }
  \item{CheckArguments}{
    Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  The originality of species is their frequency resulting in the maximum value of generalized entropy of order \eqn{q}.
  Originality was introduced by Pavoine et al. (2005) for Rao quadratic entropy (\eqn{q=2}). It can be calculated exactly by \code{\link{originality}}. For other orders of diversity, numeric optimization is required. 
}
\value{
  A named vector of frequencies.  
}
\references{
  Pavoine, S., Ollier, S. and Dufour, A.-B. (2005). Is the originality of a species measurable? \emph{Ecology Letters} 8: 579-586.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{originality}}
}
\examples{
  data(Paracou618)
  # Calculate Originality of a species according to a simple functional tree, q=1
  (Frequencies <- Originality.Species(EightSpTree, 1))
  # Plot the result
  if (require(ade4, quietly = TRUE)) {
    dotchart.phylog(EightSpTree, Frequencies)
  }
}