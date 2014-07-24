\name{IDRlsi}
\alias{IDRlsi}
\title{IDR likelihood support interval.}
\usage{
  IDRlsi(y, alpha = 0.05, k = 8, use.alpha = FALSE,
    pf = TRUE, converge = 1e-08, rnd = 3, start = NULL,
    trace.it = FALSE, iter.max = 24)
}
\arguments{
  \item{y}{Data vector c(y1, n1, y2, n2) where y are the
  positives, n are the total, and group 1 is compared to
  group 2.}

  \item{k}{Likelihood ratio criterion.}

  \item{alpha}{Complement of the confidence level.}

  \item{use.alpha}{Base choice of k on its relationship to
  alpha?}

  \item{pf}{Estimate \emph{IDR} or its complement
  \emph{PF}?}

  \item{trace.it}{Verbose tracking of the iterations?}

  \item{iter.max}{Maximum number of iterations}

  \item{converge}{Convergence criterion}

  \item{rnd}{Number of digits for rounding. Affects display
  only, not estimates.}

  \item{start}{describe here.}
}
\value{
  A \code{\link{rrsi}} object with the following elements.
  \item{estimate}{vector with point and interval estimate}
  \item{estimator}{either \emph{PF} or \emph{IDR}}
  \item{y}{data vector} \item{k}{Likelihood ratio
  criterion} \item{rnd}{how many digits to round the
  display} \item{alpha}{complement of confidence level}
}
\description{
  Estimates likelihood support interval for the incidence
  density ratio or prevented fraction based on it.
}
\details{
  Estimates likelihood support interval for the incidence
  density ratio based on orthogonal factoring of
  reparameterized likelihood. The incidence density is the
  number of cases per subject-time; its distribution is
  assumed Poisson. \cr \cr Likelihood support intervals are
  usually formed based on the desired likelihood ratio,
  often 1/8 or 1/32. Under some conditions the log
  likelihood ratio may follow the chi square distribution.
  If so, then \eqn{\alpha=1-F(2log(k),1)}, where \eqn{F} is
  a chi-square CDF. \code{RRsc()} will make the conversion
  from \eqn{\alpha} to \emph{k} if \code{use.alpha = TRUE}.
  \cr \cr The data may also be a matrix. In that case
  \code{y} would be entered as \code{matrix(c(y1, n1 - y1,
  y2, n2 - y2), 2, 2, byrow = TRUE)}.
}
\note{
  Level tested: Low.
}
\examples{
IDRlsi(c(26, 204, 10, 205), pf = FALSE)

# 1/8 likelihood support interval for IDR

# corresponds to 95.858\% confidence
#   (under certain assumptions)

# IDR
#  IDR   LL   UL
# 2.61 1.26 5.88
}
\author{
  David Siev \email{david.siev@aphis.usda.gov}
}
\references{
  Royall R. \emph{Statistical Evidence: A Likelihood
  Paradigm}. Chapman & Hall, Boca Raton, 1997. Section 7.2.
}
\seealso{
  \code{\link{IDRsc}}
}

