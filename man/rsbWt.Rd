\name{rsbWt}
\alias{rsbWt}
\title{Rao-Scott weighting.}
\usage{
  rsbWt(fit = NULL, subset.factor = NULL, fit.only = T)
}
\arguments{
  \item{fit}{A \code{\link{glm}} object.}

  \item{subset.factor}{Factor for estimating the weights by
  subset.}

  \item{fit.only}{Return only the new fit? If FALSE, also
  returns the weights and phi estimates.}
}
\value{
  A list with the following elements.  \item{fit}{the new
  model fit, updated by the estimated weights}
  \item{weights}{vector of weights} \item{d}{vector of
  \eqn{{d}_{i}}{d_i} estimates}
}
\description{
  Rao-Scott weighting of clustered binomial observations.
}
\details{
  Estimates the cluster design effect \eqn{{d}_{i}}{d_i} as
  the variance inflation due to clustering by the method of
  Rao and Scott. Observations are then weighted by the
  inverse of the \eqn{{d}_{i}}{d_i}.
}
\note{
  Level tested: Moderate.
}
\examples{
birdm.fit <- glm(cbind(y,n-y)~tx-1,binomial,birdm)
RRor(rsbWt(birdm.fit))
#
# 95\% t intervals on 4 df
#
# PF
#     PF     LL     UL
#  0.479 -1.061  0.868
#
#       mu.hat    LL     UL
# txcon  0.768 0.968 0.2659
# txvac  0.400 0.848 0.0737
#
# See the package vignette for more examples
}
\author{
  David Siev \email{david.siev@aphis.usda.gov}
}
\references{
  Rao JNK, Scott AJ, 1992. A simple method for the analysis
  of clustered binary data. \emph{Biometrics} 48:577-585.
}
\seealso{
  \code{\link{RRor}, \link{rsb}}. See the package vignette
  for more examples.
}

