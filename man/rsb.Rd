\name{rsb}
\alias{rsb}
\title{Rao-Scott weights.}
\usage{
  rsb(y, n, id = NULL)
}
\arguments{
  \item{y}{Number positive.}

  \item{n}{Total number.}

  \item{id}{Factor for estimating the weights by subset.}
}
\value{
  A list with the following elements.  \item{w}{vector of
  weights} \item{d}{vector of \eqn{{d}_{i}}{d_i} estimates}
}
\description{
  Rao-Scott weights.
}
\details{
  Estimates the cluster design effect \eqn{{d}_{i}}{d_i} as
  the variance inflation due to clustering by the method of
  Rao and Scott. \code{rsb} estimates the
  \eqn{{d}_{i}}{d_i} for use by \code{rsbWt} or other
  functions.
}
\note{
  Level tested: Moderate.
}
\examples{
# Weil's rat data (Table 1 of Rao and Scott)
rsb(rat$y, rat$n, rat$group)$d
#  control  treated
# 1.232495 3.952861
}
\author{
  David Siev \email{david.siev@aphis.usda.gov}
}
\references{
  Rao JNK, Scott AJ, 1992. A simple method for the analysis
  of clustered binary data. \emph{Biometrics} 48:577-585.
}
\seealso{
  \code{\link{rsbWt}}. See the package vignette for more
  examples.
}

