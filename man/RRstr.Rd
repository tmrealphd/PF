\name{RRstr}
\alias{RRstr}
\title{Gart-Nam method, CI for common RR over strata or clusters.}
\usage{
  RRstr(formula = NULL, data = NULL, compare = c("b", "a"),
    Y, alpha = 0.05, pf = TRUE, trace.it = FALSE,
    iter.max = 24, converge = 1e-06, rnd = 3,
    multiplier = 0.7, divider = 1.1)
}
\arguments{
  \item{formula}{Formula of the form \code{cbind(y, n) ~ x
  + cluster(w)}, where y is the number positive, n is the
  group size, x is a factor with two levels of treatment,
  and w is a factor indicating the clusters.}

  \item{data}{data.frame containing variables of formula}

  \item{compare}{Text vector stating the factor levels:
  compare[1] is the control or reference group to which
  compare[2] is compared}

  \item{Y}{Matrix of data. Each row is a stratum or
  cluster. The columns are y2, n2, y1, n1. If data entered
  by formula and dataframe, Y is generated automatically.}

  \item{pf}{Estimate \emph{RR} or its complement
  \emph{PF}?}

  \item{alpha}{Size of the homogeneity test and complement
  of the confidence level.}

  \item{trace.it}{verbose tracking of the iterations?}

  \item{iter.max}{Maximum number of iterations}

  \item{converge}{Convergence criterion}

  \item{rnd}{Number of digits for rounding. Affects display
  only, not estimates.}

  \item{multiplier}{internal control parameter for
  algorithm}

  \item{divider}{internal control parameter for algorithm}
}
\value{
  A \code{\link{rrstr}} object with the following fields:
  \item{estimate}{matrix of point and interval estimates -
  starting value, MLE, and skewness corrected}
  \item{hom}{list of homogeneity statistic, p-value, and
  degrees of freedom} \item{estimator}{either \code{"PF"}
  or \code{"RR"}} \item{y}{Y matrix of the data}
  \item{compare}{groups compared} \item{rnd}{how many
  digits to round the display} \item{alpha}{size of test;
  complement of confidence level}
}
\description{
  Estimates confidence intervals for the risk ratio or
  prevented fraction from clustered or stratified data.
}
\details{
  Uses the DUD algorithm to estimate confidence intervals
  by the method of Gart.
}
\note{
  See vignette \emph{Examples for Stratified Designs} for
  more examples (enter \code{?PF}).

  Level tested: High. \cr \cr Call to this function may be
  one of two formats: (1) specify \code{data} and
  \code{formula} or (2) as a matrix \code{Y} \cr \cr
  \code{RRstr(formula, data, compare = c('b','a'), pf =
  TRUE, alpha = 0.05, trace.it = FALSE, iter.max = 24,
  converge = 1e-6, rnd = 3, multiplier = 0.7, divider =
  1.1)} \cr \cr \code{RRstr(Y, compare = c('b','a'), pf =
  TRUE, alpha = 0.05, trace.it = FALSE, iter.max = 24,
  converge = 1e-6, rnd = 3, multiplier = 0.7, divider =
  1.1)}
}
\examples{
## Table 1 from Gart (1985)
##  as data frame
RRstr(cbind(y,n) ~ tx + cluster(clus), Table6 , pf = FALSE)

# Test of homogeneity across clusters

# stat     0.954
# df       3
# p        0.812

#  RR estimates

# 		        RR   LL   UL
# starting   2.66 1.37 5.18
# mle        2.65 1.39 5.03
# skew corr  2.65 1.31 5.08

## or as matrix
RRstr(Y = table6, pf = FALSE)
}
\author{
  David Siev \email{david.siev@aphis.usda.gov}
}
\references{
  Gart JJ, 1985. Approximate tests and interval estimation
  of the common relative risk in the combination of \eqn{2
  x 2} tables. \emph{Biometrika} 72:673-677. \cr Gart JJ,
  Nam J, 1988. Approximate interval estimation of the ratio
  of binomial parameters: a review and corrections for
  skewness. \emph{Biometrics} 44:323-338. \cr Ralston ML,
  Jennrich RI, 1978. DUD, A Derivative-Free Algorithm for
  Nonlinear Least Squares. \emph{Technometrics} 20:7-14.
}
\seealso{
  \code{\link{rrstr}}
}

