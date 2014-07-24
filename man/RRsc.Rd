\name{RRsc}
\alias{RRsc}
\title{RR score based asymptotic CI.}
\usage{
  RRsc(y = NULL, alpha = 0.05, pf = TRUE, trace.it = FALSE,
    iter.max = 18, converge = 1e-06, rnd = 3)
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
}
\value{
  A \code{\link{rrsc}} object with the following fields.
  \item{estimate}{matrix of point and interval estimates -
  see details} \item{estimator}{either \code{"PF"} or
  \code{"RR"}} \item{y}{data vector} \item{rnd}{how many
  digits to round the display} \item{alpha}{complement of
  confidence level}
}
\description{
  Estimates confidence intervals for the risk ratio or
  prevented fraction based on the score statistic.
}
\details{
  Estimates are returned for three estimators based on the
  score statistic. The score method was introduced by
  Koopman (1984). Gart and Nam's modification (1988)
  includes a skewness correction. The method of Miettinen
  and Nurminen (1985) is a version made slightly more
  conservative than Koopman's by including a factor of
  \code{(N-1)/N}. The starting estimate for the DUD
  algorithm is obtained by the modified Katz method (log
  method with 0.5 added to each cell). Both forms of the
  Katz estimate may be retrieved from the returned object
  using \code{RRsc()$estimate}. \cr \cr The data may also
  be a matrix. In that case \code{y} would be entered as
  \code{matrix(c(y1, n1-y1, y2, n2-y2), 2, 2, byrow =
  TRUE)}.
}
\note{
  Level tested: High.
}
\examples{
RRsc(c(4, 24, 12, 28))

# PF
# 95\% interval estimates

# 				PF     LL    UL
# MN method    0.611 0.0251 0.857
# score method 0.611 0.0328 0.855
# skew corr    0.611 0.0380 0.876
}
\author{
  David Siev \email{david.siev@aphis.usda.gov}
}
\references{
  Gart JJ, Nam J, 1988. Approximate interval estimation of
  the ratio of binomial parameters: a review and
  corrections for skewness. \emph{Biometrics} 44:323-338.
  \cr Koopman PAR, 1984. Confidence intervals for the ratio
  of two binomial proportions. \emph{Biometrics}
  40:513-517. \cr Miettinen O, Nurminen M, 1985.
  Comparative analysis of two rates. \emph{Statistics in
  Medicine} 4:213-226. \cr Ralston ML, Jennrich RI, 1978.
  DUD, A Derivative-Free Algorithm for Nonlinear Least
  Squares. \emph{Technometrics} 20:7-14.
}
\seealso{
  \code{\link{rrsc}}
}

