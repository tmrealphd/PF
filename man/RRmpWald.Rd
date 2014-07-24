\name{RRmpWald}
\alias{RRmpWald}
\title{Wald confidence intervals for RR from matched pairs}
\usage{
  RRmpWald(formula = NULL, data = NULL,
    compare = c("con", "vac"), affected = 1, x,
    alpha = 0.05, pf = TRUE, tdist = TRUE, df = NULL,
    rnd = 3)
}
\arguments{
  \item{formula}{Formula of the form \code{y ~ x +
  cluster(w)}, where y is the indicator for an individual's
  positive response, x is a factor with two levels of
  treatment, and w identifies the pairs.}

  \item{data}{\code{data.frame} containing variables in
  formula}

  \item{compare}{Text vector stating the factor levels:
  compare[1] is the control or reference group to which
  compare[2] is compared}

  \item{affected}{Indicator for positive response}

  \item{x}{Alternative data input. Instead of formula and
  data frame, data may be input as vector or 2x2 table}

  \item{alpha}{Complement of the confidence level}

  \item{pf}{Estimate \emph{RR} or its complement
  \emph{PF}?}

  \item{tdist}{Use t distribution?}

  \item{df}{Degrees of freedom. When NULL, the function
  will default to \code{df = N - 2}, where N is the total
  number of pairs.}

  \item{rnd}{Number of digits for rounding. Affects display
  only, not estimates.}
}
\value{
  A \code{\link{rrmp}} object with the following fields:
  \item{estimate}{vector of point and interval estimates -
  see details} \item{estimator}{either \code{"PF"} or
  \code{"RR"}} \item{compare}{text vector, same as input}
  \item{alpha}{complement of confidence level}
  \item{rnd}{how many digits to round the display}
  \item{xtable}{data arrayed in 2x2 matrix}
  \item{freqvec}{data arrayed a 4-vector}
  \item{multvec}{data frame showing the multinomial
  representation of the data}
}
\description{
  Estimates confidence intervals for the risk ratio or
  prevented fraction from matched pairs.
}
\details{
  Estimates confidence intervals for the risk ratio or
  prevented fraction from matched pairs. The response is
  the tetranomial vector [11, 12, 21, 22], where the first
  index is the row and the the second index is the column
  when displayed as a 2x2 table. Wald type confidence
  intervals are found by applying the delta method to the
  multinomial variance. This method fails when there are no
  responders in one of the treatment groups. \cr \cr
  Alternative forms of data entry are illustrated by the
  output, say \code{Y}, where \code{c(Y$xtable) = Y$freqvec
  = Y$multvec$Freq}. \cr \cr If RR = 0 (PF = 1), the
  function will return degenerate interval.
}
\note{
  Level tested: Low. \cr\cr Experimental functions for
  estimating profile likelihood intervals are in the
  CVBmisc package. \cr \cr Call to this function may be one
  of two formats: (1) specify \code{data} and
  \code{formula} or (2) as a vector or 2 x 2 table \code{x}
  \cr \cr \code{RRmpWald(formula, data, compare = c('con',
  'vac'), affected = 1, alpha = 0.05, pf = TRUE, tdist =
  TRUE, df = NULL, rnd = 3)} \cr \cr \code{RRmpWald(x,
  compare = c('con', 'vac'), affected = 1, alpha = 0,05, pf
  = TRUE, tdist = TRUE, df = NULL, rnd = 3)}
}
\examples{
RRmpWald(pos ~ tx + cluster(cage), New, compare = c('con', 'vac'))

# PF
# 95\% interval estimates
#
#   PF    LL    UL
# 0.550 0.183 0.752


RRmpWald(x = c(2, 9, 1, 6))

# PF
# 95\% interval estimates

#   PF    LL    UL
# 0.727 0.124 0.915
}
\author{
  David Siev \email{david.siev@aphis.usda.gov}
}

