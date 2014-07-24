\name{tauWt}
\alias{tauWt}
\title{Binomial dispersion: intra-cluster correlation parameter.}
\usage{
  tauWt(fit, subset.factor = NULL, fit.only = TRUE,
    iter.max = 12, converge = 1e-06, trace.it = FALSE)
}
\arguments{
  \item{fit}{A \code{\link{glm}} object.}

  \item{subset.factor}{Factor for estimating tau by
  subset.}

  \item{fit.only}{Return only the final fit?  If FALSE,
  also returns the weights and tau estimates.}

  \item{iter.max}{Maximum number of iterations.}

  \item{converge}{Convergence criterion: difference between
  model degrees of freedom and Pearson's chi-square.
  Default 1e-6.}

  \item{trace.it}{Display print statments indicating
  progress}
}
\value{
  A list with the following elements.  \item{fit}{the new
  model fit, updated by the estimated weights}
  \item{weights}{vector of weights} \item{phi}{vector of
  phi estimates}
}
\description{
  MME estimates of binomial dispersion parameter tau
  (intra-cluster correlation).
}
\details{
  Estimates binomial dispersion parameter \eqn{\tau} by the
  method of moments. Iteratively refits the model by the
  Williams procedure, weighting the observations by
  \eqn{1/\phi_{ij}}{1/\phi_ij}, where \eqn{\phi_{ij}=1+\tau
  _j(n_{ij}-1)}{\phi_ij=1+\tau_j(n_ij - 1)}, \eqn{j}
  indexes the subsets, and \eqn{i} indexes the
  observations.
}
\note{
  Level tested: Moderate. (Extensive in S+, not as much in
  R. They have different \code{glm()} functions.)
}
\examples{
birdm.fit <- glm(cbind(y,n-y)~tx-1, binomial, birdm)
RRor(tauWt(birdm.fit))

# 95\% t intervals on 4 df
#
# PF
#     PF     LL     UL
#  0.489 -0.578  0.835
#
#       mu.hat    LL    UL
# txcon  0.737 0.944 0.320
# txvac  0.376 0.758 0.104
#
# See the package vignette for more examples
}
\author{
  David Siev \email{david.siev@aphis.usda.gov}
}
\references{
  Williams DA, 1982. Extra-binomial variation in logistic
  linear models. \emph{Applied Statistics} 31:144-148. \cr
  Wedderburn RWM, 1974. Quasi-likelihood functions,
  generalized linear models, and the Gauss-Newton method.
  \emph{Biometrika} 61:439-447.
}
\seealso{
  \code{\link{phiWt}}, \code{\link{RRor}}, see package
  vignette for more examples
}

