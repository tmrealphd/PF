\name{RRtosst}
\alias{RRtosst}
\title{RR exact CI, TOSST method.}
\usage{
  RRtosst(y, alpha = 0.05, pf = TRUE, stepstart = 0.1,
    iter.max = 36, converge = 1e-06, rnd = 3,
    trace.it = FALSE, nuisance.points = 120, gamma = 1e-06)
}
\arguments{
  \item{y}{Data vector c(y1, n1, y2, n2) where y are the
  positives, n are the total, and group 1 is compared to
  group 2.}

  \item{alpha}{Complement of the confidence level.}

  \item{pf}{Estimate \emph{RR} or its complement
  \emph{PF}?}

  \item{trace.it}{Verbose tracking of the iterations?}

  \item{iter.max}{Maximum number of iterations}

  \item{converge}{Convergence criterion}

  \item{rnd}{Number of digits for rounding. Affects display
  only, not estimates.}

  \item{stepstart}{starting interval for step search}

  \item{nuisance.points}{number of points over which to
  evaluate nuisance parameter}

  \item{gamma}{parameter for Berger-Boos correction
  (restricts range of nuisance parameter evaluation)}
}
\value{
  A \code{\link{rr1}} object with the following fields.
  \item{estimate}{vector with point and interval estimate}
  \item{estimator}{either \code{"PF"} or \code{"RR"}}
  \item{y}{data vector} \item{rnd}{how many digits to round
  the display} \item{alpha}{complement of confidence level}
}
\description{
  Estimates confidence interval for the risk ratio or
  prevented fraction; exact method based on the score
  statistic (inverts two one-sided tests).
}
\details{
  Estimates confidence intervals based on the score
  statistic that are 'exact' in the sense of accounting for
  discreteness. Inverts two one-sided score tests. The
  score statistic is used to select tail area tables, and
  the binomial probability is estimated over the tail area
  by taking the maximum over the nuisance parameter.
  Algorithm is a simple step search. \cr \cr The data may
  also be a matrix. In that case \code{y} would be entered
  as \code{matrix(c(y1, n1-y1, y2, n2-y2), 2, 2, byrow =
  TRUE)}.
}
\note{
  Level tested: Moderate.
}
\examples{
\dontrun{RRtosst(c(4, 24, 12, 28))

# PF
# 95\% interval estimates

#    PF    LL    UL
# 0.611 0.012 0.902 }
}
\author{
  David Siev \email{david.siev@aphis.usda.gov}
}
\references{
  Koopman PAR, 1984. Confidence intervals for the ratio of
  two binomial proportions. \emph{Biometrics} 40:513-517.
  \cr Agresti A, Min Y, 2001.  On small-sample confidence
  intervals for parameters in discrete distribution.
  \emph{Biometrics} 57: 963-971. \cr Berger RL, Boos DD,
  1994. P values maximized over a confidence set for the
  nuisance parameter. \emph{Journal of the American
  Statistical Association} 89:214-220.
}
\seealso{
  \code{\link{RRotsst}, \link{rr1}}
}

