\name{RRor}
\alias{RRor}
\title{RR estimate from logistic regression.}
\usage{
  RRor(fit = NULL, beta.hat = NULL, var.beta.hat = NULL,
    degf = NULL, which = c(1, 2), pf = T, norm = F,
    alpha = 0.05, rnd = 3)
}
\arguments{
  \item{fit}{A \code{\link{glm}} object.}

  \item{beta.hat}{Parameters estimates from a logistic
  regression with no intercept.}

  \item{var.beta.hat}{Variance-covariance matrix from a
  logistic regression with no intercept.}

  \item{degf}{Degrees of freedom.}

  \item{which}{Numeric vector indicating which parameters
  to compare, so that \emph{RR} = compare[2]/compare[1]}

  \item{pf}{Estimate \emph{RR} or its complement
  \emph{PF}?}

  \item{norm}{Estimate confidence interval using quantiles
  of Guassian rather than t distribution quantiles?}

  \item{alpha}{Complement of the confidence level.}

  \item{rnd}{Number of digits for rounding. Affects display
  only, not estimates.}
}
\value{
  A \code{\link{rror}} object with the following fields.
  \item{estimate}{vector with point and interval estimate}
  \item{estimator}{either \emph{PF} or \emph{RR}}
  \item{mu}{matrix with rows giving probability estimates
  for each of the groups} \item{rnd}{how many digits to
  round the display} \item{alpha}{complement of confidence
  level} \item{norm}{logical indicating Gaussian or t
  interval} \item{degf}{degrees of freedom}
}
\description{
  Model based interval estimate of the risk ratio or
  prevented fraction from a logistic regression model.
}
\details{
  Estimates confidence intervals using the delta method on
  parameters from a generalized linear model with logit
  link. \cr \cr
  \eqn{RR={{{\mu}}_{2}}/{{{\mu}}_{1}}}{RR=muhat_2/muhat_1},
  where \eqn{{\mu}_{i}}{muhat_i} are the estimated
  probabilities from the model.
}
\note{
  Level tested: Moderate. (Extensive in S+, not as much in
  R. They have different \code{glm()} functions.) \cr \cr
  Call to this function may be one of two formats: (1)
  specify \code{fit} or (2) \code{beta.hat},
  \code{var.beta.hat}, \code{degf} \cr \cr \code{RRor(fit,
  degf = NULL, pf = TRUE, alpha = 0.05, which = c(1, 2),
  norm = TRUE, rnd = 3)} \cr \cr \code{RRor(beta.hat,
  var.beta.hat, degf, pf = TRUE, alpha = 0.05, which = c(1,
  2), norm = TRUE, rnd = 3)}
}
\examples{
bird.fit <- glm(cbind(y, n - y) ~ tx - 1, binomial, bird)
RRor(tauWt(bird.fit))

# 95\% t intervals on 4 df
#
# PF
#     PF     LL     UL
#  0.500 -0.583  0.842
#
#       mu.hat    LL     UL
# txcon  0.733 0.943 0.3121
# txvac  0.367 0.752 0.0997
#
# See the package vignette for more examples
}
\author{
  David Siev \email{david.siev@aphis.usda.gov}
}
\references{
  No references yet.
}
\seealso{
  \code{\link{rror}, \link{phiWt}, \link{tauWt}}. See the
  package vignette for more examples.
}

