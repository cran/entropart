\name{expq}
\alias{expq}
\title{
  Exponential of order q
}
\description{
  Calculates the deformed exponential of order \eqn{q}.
}
\usage{
expq(x, q)
}
\arguments{
  \item{x}{
    A numeric vector.
  }
  \item{q}{
    A number.
  }
}
\details{
  The deformed exponential is defined as \eqn{(x(1-q)+1)^{\frac{1}{(1-q)}}}.
  
  For \eqn{q>1}, \eqn{\ln_q{(+\infty)}=\frac{1}{(q-1)}} so \eqn{\exp_q{(x)}} is not defined for \eqn{x>\frac{1}{(q-1)}}.
}
\value{
  A vector of the same length as \code{x} containing the transformed values.
}
\references{
  Marcon, E., Scotti, I., Herault, B., Rossi, V. and Lang, G. (2014). Generalization of the partitioning of Shannon diversity. \emph{PLOS One} 9(3): e90289.
  
  Tsallis, C. (1994). What are the numbers that experiments provide? \emph{Quimica Nova} 17(6): 468-471.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{expq}}
}
\examples{
  curve(exp(x), -5, 0, lty=3)  
  curve(expq(x, 2), -5, 0, lty=2, add=TRUE)
  curve(expq(x, 3), -5, 0, lty=1, add=TRUE)
  legend("topleft", legend = c("exp(x)", "exp2(x)", "exp3(x)"), lty = c(1, 2, 3), inset=0.02)
}