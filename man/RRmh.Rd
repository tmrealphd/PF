\name{RRmh}
\alias{RRmh}
\title{Mantel-Haenszel method, CI for common RR over strata or clusters with sparse data.}
\usage{
  RRmh(formula = NULL, data = NULL, compare = c("b", "a"),
    Y, alpha = 0.05, pf = TRUE, rnd = 3)
}
\arguments{
  \item{formula}{Formula of the form \code{cbind(y, n) ~ x
  + cluster(w)}, where \code{y} is the number positive,
  \code{n} is the group size, \code{x} is a factor with two
  levels of treatment, and \code{w} is a factor indicating
  the clusters.}

  \item{data}{\code{data.frame} containing variables for
  formula}

  \item{compare}{Text vector stating the factor levels:
  \code{compare[1]} is the control or reference group to
  which \code{compare[2]} is compared}

  \item{Y}{Matrix of data, \eqn{K \times 4}{K x 4}. Each
  row is a stratum or cluster. The columns are \eqn{y2, n2,
  y1, n1}, where the y's are the number of positive in each
  group, and the n is the total in each group. If data
  entered by formula and dataframe, \code{Y} is generated
  automatically.}

  \item{pf}{Estimate \emph{RR} or its complement
  \emph{PF}?}

  \item{alpha}{Complement of the confidence level.}

  \item{rnd}{Number of digits for rounding. Affects display
  only, not estimates.}
}
\value{
  An object of class \code{\link{rr1}} with the following
  fields.  \item{estimate}{vector of point and interval
  estimates: point estimate, lower confidence limit, upper
  confidence limit} \item{estimator}{either \code{"PF"} or
  \code{"RR"}} \item{y}{y matrix of the data}
  \item{compare}{groups compared} \item{rnd}{how many
  digits to round the display} \item{alpha}{complement of
  confidence level}
}
\description{
  Estimates confidence intervals for the risk ratio or
  prevented fraction from clustered or stratified data,
  using a Mantel-Haenszel estimator for sparse data.
}
\details{
  Based on the Mantel-Haenszel (1959) procedure for sparse
  data developed by Greenland and Robins (1985). The
  confidence limits are based on asymptotic normality of
  the log(risk ratio).  Agresti and Hartzel (2000) favor
  this procedure for small, sparse data sets, but they warn
  that it is less efficient than maximum likelihood for
  large data sets.
}
\note{
  If either all y1's or all y2's are zero, a division by
  zero may occur, and a NaN returned for some values. \cr
  \cr Level tested: Low. \cr\cr See vignette \emph{Examples
  for Stratified Designs} for more examples (enter
  \code{?PF}). \cr Call to this function may be one of two
  formats: (1) specify \code{data} and \code{formula} or
  (2) as a matrix \code{Y} \cr \cr \code{RRmh(formula,
  data, compare = c('b','a'), pf = TRUE, alpha = 0.05, rnd
  = 3)} \cr \cr \code{RRmh(Y, pf = TRUE, alpha = 0.05, rnd
  = 3)}
}
\examples{
## Table 1 from Gart (1985)
##  as data frame
RRmh(cbind(y,n) ~ tx + cluster(clus), Table6 , pf = FALSE)

#  RR estimates

# RR
# 95\% interval estimates
#
#   RR   LL   UL
# 2.67 1.37 5.23
#

## or as matrix
RRmh(Y = table6, pf = FALSE)
}
\author{
  Christopher Tong
  \email{Christopher.H.Tong@aphis.usda.gov}
}
\references{
  Mantel N, Haenszel W, 1959.  Statistical aspects of the
  analysis of data from retrospective studies of disease.
  \emph{Journal of the National Cancer Institute} 22:
  719-748. \cr \cr Greenland S, Robins JM, 1985.
  Estimation of a common effect parameter from sparse
  follow-up data.  \emph{Biometrics} 41: 55-68.  Errata,
  45: 1323-1324. \cr \cr Agresti A, Hartzel J, 2000.
  Strategies for comparing treatments on a binary response
  with multi-centre data.  \emph{Statistics in Medicine}
  19: 1115-1139. \cr Lachin JM, 2000.  \emph{Biostatistical
  Methods: The Assessment of Relative Risks} (Wiley, New
  York), Sec. 4.3.1.
}
\seealso{
  \code{\link{rr1}}
}

