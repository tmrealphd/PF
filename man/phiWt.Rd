\name{phiWt}
\alias{phiWt}
\title{Binomial dispersion parameter.}
\usage{
  phiWt(fit, subset.factor = NULL, fit.only = TRUE,
    show.warns = FALSE)
}
\arguments{
  \item{fit}{A \code{\link{glm}} object.}

  \item{subset.factor}{Factor for estimating phi by
  subset.}

  \item{fit.only}{Return only the new fit?  If FALSE, also
  returns the weights and phi estimates.}

  \item{show.warns}{Show warnings}
}
\value{
  A list with the following elements.  \item{fit}{the new
  model fit, updated by the estimated weights}
  \item{weights}{vector of weights} \item{phi}{vector of
  phi estimates}
}
\description{
  MME estimate of dispersion parameter phi.
}
\details{
  Estimates binomial dispersion parameter \eqn{\phi} by the
  method of moments. Refits the model, weighting the
  observations by \eqn{1/\phi}. Uses \code{quasibinomial}
  family in \code{glm()}.
}
\note{
  Level tested: Moderate. (Extensive in S+, not as much in
  R. They have different \code{glm()} functions.)
}
\examples{
birdm.fit <- glm(cbind(y, n - y)~tx-1, binomial, birdm)
RRor(phiWt(birdm.fit))
#
# 95\% t intervals on 4 df
#
# PF
#     PF     LL     UL
#  0.479 -0.537  0.823
#
#       mu.hat   LL    UL
# txcon  0.768 0.95 0.367
# txvac  0.400 0.78 0.111
#
# See the package vignette for more examples
}
\author{
  David Siev \email{david.siev@aphis.usda.gov}
}
\references{
  Wedderburn RWM, 1974. Quasi-likelihood functions,
  generalized linear models, and the Gauss-Newton method.
  \emph{Biometrika} 61:439-447.
}
\seealso{
  \code{\link{tauWt}}, \code{\link{RRor}}. See the package
  vignette for more examples.
}

