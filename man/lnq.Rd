\name{lnq}
\alias{lnq}
\title{
  Logarithm of order q
}
\description{
  Calculates the deformed logarithm of order \eqn{q}.
}
\usage{
lnq(x, q)
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
  The deformed logarithm is defined as \eqn{\ln_q{x}=\frac{(x^{(1-q)}-1)}{(1-q)}}.
  
  The shape of the deformed logarithm is similar to that of the regular one. \eqn{\ln_1{x}=\log{x}}.
  
  For \eqn{q>1}, \eqn{\ln_q{(+\infty)}=\frac{1}{(q-1)}}.
}
\value{
  A vector of the same length as \code{x} containing the transformed values.  
}
\references{
  Tsallis, C. (1994). What are the numbers that experiments provide? \emph{Quimica Nova} 17(6): 468-471.
}
\author{
  Eric Marcon <Eric.Marcon@ecofog.gf>
}
\seealso{
  \code{\link{expq}}
}
\examples{
  curve(log(x), 0, 1, lty=1)
  curve(lnq(x, 2), 0, 1, lty=2, add=TRUE)
  curve(lnq(x, 3), 0, 1, lty=3, add=TRUE)  
  legend("topleft", legend = c("log(x)", "ln2(x)", "ln3(x)"), lty = c(1, 2, 3), inset=0.02)
}