\name{CommunityProfile}
\alias{CommunityProfile}
\alias{as.CommunityProfile}
\alias{is.CommunityProfile}
\alias{plot.CommunityProfile}
\alias{autoplot.CommunityProfile}
\alias{CEnvelope}
\title{
  Diversity or Entropy Profile of a community
}
\description{
  Calculates the diversity or entropy profile of a community, applying a community function to a vector of orders.
}
\usage{
CommunityProfile(FUN, NorP, q.seq = seq(0, 2, 0.1), 
    NumberOfSimulations = 0, Alpha = 0.05, BootstrapMethod = "Chao2015", 
    size = 1, \dots, ShowProgressBar = TRUE, CheckArguments = TRUE)
as.CommunityProfile(x, y, low = NULL, high = NULL, mid = NULL)
is.CommunityProfile(x)
\method{plot}{CommunityProfile}(x, \dots, main = NULL, 
          xlab = "Order of Diversity", ylab = "Diversity", ylim = NULL,
          LineWidth = 2, ShadeColor = "grey75", BorderColor = "red")
\method{autoplot}{CommunityProfile}(object, \dots, main = NULL, 
          xlab = "Order of Diversity", ylab = "Diversity",
          ShadeColor = "grey75", alpha = 0.3, BorderColor = "red",
          col = ggplot2::GeomLine$default_aes$colour,
          lty = ggplot2::GeomLine$default_aes$linetype,
          lwd = ggplot2::GeomLine$default_aes$size)
CEnvelope(Profile, LineWidth = 2, ShadeColor = "grey75", BorderColor = "red", \dots)
}
\arguments{
  \item{FUN}{
  The function to be applied to each value of \eqn{q.seq}. Any function accepting a numeric vector (or a two-column matrix) and a number as first two arguments and an argument named \code{CheckArguments} is acceptable (other arguments of the functions are passed by \code{\dots}). See *Details* for useful entropy and diversity functions and *Examples* for an ad-hoc one.
  }
  \item{NorP}{
  A numeric vector. Contains either abundances or probabilities.
  }
  \item{q.seq}{
  A numeric vector: the sequence of diversity orders to address. Default is from 0 to 2.
  }
  \item{NumberOfSimulations}{
  The number of simulations to run, 0 by default.
  }
  \item{Alpha}{
  The risk level, 5\% by default.
  }
  \item{BootstrapMethod}{
  The method used to obtain the probabilities to generate bootstrapped communities from observed abundances. See \code{\link{rCommunity}}.
  }
  \item{size}{
  The size of simulated communities used to compute the bootstrap confidence envelope. 1 (default) means that the actual size must be used.
  }
  \item{object}{
  An object.
  }
  \item{x}{
  An object to be tested or plotted or the vector of orders of community profiles in \code{as.CommunityProfile}.
  }
  \item{y}{
  Entropy or diversity values of each order, corresponding to \code{x} values.
  }
  \item{low}{
  Entropy or diversity lower bound of the confidence envelope, corresponding to \code{x} values.
  }
  \item{high}{
  Entropy or diversity higher bound of the confidence envelope, corresponding to \code{x} values.
  }
  \item{mid}{
  Entropy or diversity center value (usually the mean) of the confidence envelope, corresponding to \code{x} values.
  }
  \item{Profile}{
  An \code{CommunityProfile} to be plotted.
  }
  \item{\dots}{
  Additional arguments to be passed to \code{FUN} in \code{CommunityProfile}, to \code{plot} in \code{plot.CommunityProfile} or to \code{lines} in \code{CEnvelope}.
  }
  \item{main}{
  The main title of the plot.
  }
  \item{xlab}{
  The x axis label of the plots.
  }
  \item{ylab}{
  The y axis label of the plot.
  }
  \item{ylim}{
  The interval of y values plotted.
  }
  \item{LineWidth}{
  The width of the line that represents the actual profile.
  }
  \item{ShadeColor}{
  The color of the shaded confidence envelope.
  }
  \item{BorderColor}{
  The color of the bounds of the confidence envelope.
  }
  \item{alpha}{
  Opacity of the confidence enveloppe, between 0 and 1.
  }
  \item{col}{
  The color of the geom objects. See "Color Specification" in \code{\link{par}}.
  }
  \item{lty}{
  The type of the lines. See \code{\link{lines}}.
  }
  \item{lwd}{
  The width of the lines. See \code{\link{lines}}.
  }
  \item{ShowProgressBar}{
  If \code{TRUE} (default), a progress bar is shown.
  }
  \item{CheckArguments}{
  Logical; if \code{TRUE}, the function arguments are verified. Should be set to \code{FALSE} to save time when the arguments have been checked elsewhere.
  }
}
\details{
  The function \code{CommunityProfile} is used to calculate diversity or entropy profiles based on community functions such as \code{\link{Tsallis}} or \code{\link{ChaoPD}}. 
  The first two arguments of the function must be a probability or abundance vector and a number (\eqn{q}).
  Additional arguments cannot be checked. 
  Unexpected results may be returned if \code{FUN} is not used properly.
  
  If \code{NumberOfSimulations} is greater than 0, a bootstrap confidence interval is produced by simulating communities with \code{\link{rCommunity}} and calculating their profiles. 
  The size of those communities may be that of the actual community or specified by \code{size}.
  Simulating communities implies a downward bias in the estimation: rare species of the actual community may have abundance zero in simulated communities.
  Simulated diversity values are recentered if `size = 1` so that their mean is that of the actual community.
  Else, it is assumed that the bias is of interest and must not be corrected.
  
  \code{CommunityProfile} objects can be plotted. 
  They can also be added to the current plot by \code{CEnvelope}.
}
\value{
  A \code{CommunityProfile}, which is a list:
  \item{x}{The order \eqn{q} values}
  \item{y}{The entropy or diversity values returned by \code{FUN}}
  \item{low}{The lower bound of the confidence interval}
  \item{high}{The upper bound of the confidence interval}
}
\author{
  Eric Marcon <Eric.Marcon@agroparistech.fr>, Bruno Herault <Bruno.Herault@cirad.fr>
}
\examples{
# Load Paracou data (number of trees per species in two 1-ha plot of a tropical forest)
data(Paracou618)
# Plot diversity estimated without bias correction
plot(CommunityProfile(Diversity, Paracou618.MC$Ps, seq(0, 2, 0.2)), 
lty=3, ylim=c(50, 350))
# Estimate diversity, with a condidence envelope 
# (only 10 simulations to save time, should be 1000)
Profile <- CommunityProfile(Diversity, as.AbdVector(Paracou618.MC$Ns), 
		  seq(0, 2, 0.2), Correction="UnveilJ", NumberOfSimulations=10)
# Complete the plot, and add the legend
CEnvelope(Profile, main="Paracou Plots Diversity")
legend("topright", c("Bias Corrected", "Biased"), lty=c(1,3), inset=0.01)

# Advanced use with beta-diversity functions :
# Profile of the beta entropy of the first community of Paracou618.
# Observed and expected probabilities are bound into a 2-column matrix
# An intermediate function is necessary to separate them before calling TsallisBeta
# The CheckArguments is mandatory but does not need to be set: CommunityProfile() sets it to FALSE
CommunityProfile(function(PandPexp, q, CheckArguments) 
  {TsallisBeta(PandPexp[, 1], PandPexp[, 2], q)}, 
  NorP=cbind(Paracou618.MC$Psi[, 1], Paracou618.MC$Ps), q.seq=seq(0, 2, 0.2))
}
