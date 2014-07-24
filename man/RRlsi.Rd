\name{RRlsi}
\alias{RRlsi}
\title{RR likelihood support interval.}
\usage{
  RRlsi(y, alpha = 0.05, k = 8, use.alpha = FALSE,
    pf = TRUE, iter.max = 50, converge = 1e-06, rnd = 3,
    start = NULL, track = FALSE, full.track = FALSE)
}
\arguments{
  \item{y}{Data vector c(y1, n1, y2, n2) where y are the
  positives, n are the total, and group 1 is compared to
  group 2.}

  \item{k}{Likelihood ratio criterion.}

  \item{alpha}{Complement of the confidence level (see
  details).}

  \item{use.alpha}{Base choice of k on its relationship to
  alpha?}

  \item{pf}{Estimate \emph{RR} or its complement
  \emph{PF}?}

  \item{iter.max}{Maximum number of iterations}

  \item{converge}{Convergence criterion}

  \item{rnd}{Number of digits for rounding. Affects display
  onlyRR, not estimates.}

  \item{start}{Optional starting value.}

  \item{track}{Verbose tracking of the iterations?}

  \item{full.track}{Verbose tracking of the iterations?}
}
\value{
  An object of class \code{\link{rrsi}} with the following
  fields: \cr \item{estimate}{matrix of point and interval
  estimates - see details} \item{estimator}{either
  \code{"PF"} or \code{"RR"}} \item{y}{data vector}
  \item{rnd}{how many digits to round the display}
  \item{k}{likelihood ratio criterion}
  \item{alpha}{complement of confidence level}
}
\description{
  likelihood support interval for the risk ratio or
  prevented fraction by the likelihood profile.
}
\details{
  Estimates a likelihood support interval for \emph{RR} or
  \emph{PF} by the profile likelihood method using the DUD
  algorithm. \cr \cr Likelihood support intervals are
  usually formed based on the desired likelihood ratio,
  often 1/8 or 1/32. Under some conditions the log
  likelihood ratio may follow the chi square distribution.
  If so, then \eqn{\alpha=1-F(2log(k),1)}, where \eqn{F} is
  a chi-square CDF. \code{RRlsi()} will make the conversion
  from \eqn{\alpha} to \emph{k} if \code{use.alpha = TRUE}.
  \cr \cr The data may also be a matrix. In that case
  \code{y} would be entered as \code{matrix(c(y1, n1-y1,
  y2, n2-y2), 2, 2, byrow = TRUE)}.
}
\note{
  Level tested: High.
}
\examples{
RRlsi(c(4,24,12,28))

# 1/8 likelihood support interval for PF

# corresponds to 95.858\% confidence
#   (under certain assumptions)

# PF
#     PF     LL     UL
# 0.6111 0.0168 0.8859
}
\author{
  David Siev \email{david.siev@aphis.usda.gov}
}
\references{
  Royall R. \emph{Statistical Evidence: A Likelihood
  Paradigm}. Chapman & Hall, Boca Raton, 1997.  Section 7.6
  \cr Ralston ML, Jennrich RI, 1978. DUD, A Derivative-Free
  Algorithm for Nonlinear Least Squares.
  \emph{Technometrics} 20:7-14.
}

