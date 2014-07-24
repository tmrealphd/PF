#' @title IDR confidence interval.
#' @description Estimates confidence interval for the incidence density ratio or prevented fraction based on it.
#' @details The incidence density is the number of cases per subject-time; its distribution is assumed Poisson. \code{IDRsc} estimates 
#' a confidence interval for the incidence density ratio using Siev's formula based on the Poisson score statistic.
#' \eqn{
#' IDR=\widehat{IDR}\left\{ 1+\left( \frac{1}{{{y}_{1}}}+\frac{1}{{{y}_{2}}} \right)\frac{z_{\alpha /2}^{2}}{2}\ \ \pm \ \ \frac{z_{\alpha
#'  /2}^{2}}{2{{y}_{1}}{{y}_{2}}}\sqrt{{{y}_{\bullet }}\left( {{y}_{\bullet }}z_{\alpha /2}^{2}+4{{y}_{1}}{{y}_{2}} \right)} \right\}
#' }{(See the PF Package Vignette for the formula.)}
#' \cr \cr The data may also be a matrix. In that case \code{y} would be entered as \code{matrix(c(y1, n1 - y1, y2, n2 - y2), 2, 2, byrow = TRUE)}.
#' @param y Data vector c(y1, n1, y2, n2) where y are the positives, n are the total, and group 1 is compared to group 2.
#' @param alpha Complement of the confidence level.
#' @param pf Estimate \emph{IDR}, or its complement \emph{PF}? 
#' @param rnd Number of digits for rounding. Affects display only, not estimates.
#' @return A \code{\link{rr1}} object with the following elements.
#'  \item{estimate}{vector with point and interval estimate}
#'  \item{estimator}{either \emph{PF} or \emph{IDR}}
#'  \item{y}{data vector}
#'  \item{rnd}{how many digits to round the display}
#'  \item{alpha}{complement of confidence level}
#' @export
#' @references Siev D, 1994. Estimating vaccine efficacy in prospective studies. \emph{Preventive Veterinary Medicine} 20:279-296, Appendix 1.
#' \cr Graham PL, Mengersen K, Morton AP, 2003. Confidence limits for the ratio of two rates based on likelihood scores:non-iterative method 
#' \emph{Statistics in Medicine} 22:2071-2083.
#' \cr Siev D, 2004. Letter to the editor. \emph{Statistics in Medicine} 23:693. (Typographical error in formula: replace the two final minus
#' signs with subscript dots.)
#' @author David Siev \email{david.siev@@aphis.usda.gov}
#' @note Level tested: High.
#' @seealso \code{\link{IDRlsi}}
#' 
#' @examples
#' IDRsc(c(26, 204, 10, 205),pf = FALSE)
#' 
#' # IDR 
#' # 95% interval estimates
#' 
#' #  IDR   LL   UL 
#' # 2.61 1.28 5.34 

#-------------------------------
# IDRsc
#-------------------------------
IDRsc <- function(y, alpha = 0.05, pf = TRUE, rnd =3){
	# confidence intervals for idr by score method
	# Siev 1993

	# data vector
	if(is.matrix(y))
		y <- c(t(cbind(y[,1],apply(y,1,sum))))
	y1 <- y[1.]
	s1 <- y[2.]
	y2 <- y[3.]
	s2 <- y[4.]
	y.dot <- y1 + y2
	z <- qnorm(1. - alpha/2.)
	idr.hat <- (y1 * s2)/(y2 * s1)
	det <- (z * sqrt(y.dot * (y.dot * z^2. + 4. * y1 * y2)))/(2. * y1 * y2)
	det <- c( - det, det)
	ci <- idr.hat * (1. + (z^2. * (1./y1 + 1./y2))/2. + det)
	int <- c(idr.hat, ci)
	if(!pf) 
        names(int) <- c("IDR", "LL", "UL")
    else{
        int <- 1 - int[c(1,3,2)]
        names(int) <- c("PF.IDR", "LL", "UL")
        }
	return(rr1$new(estimate = int, estimator = ifelse(pf, 'PF_IDR', 'IDR'), y = as.matrix(y), rnd = rnd, alpha = alpha))
    # out <- list(estimate = int, estimator = ifelse(pf, 'PF_IDR', 'IDR'), y = y, rnd = rnd, alpha = alpha)
    # class(out) <- 'rr1'
    # return(out)
}
