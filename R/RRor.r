#' @title RR estimate from logistic regression. 
#' @description Model based interval estimate of the risk ratio or prevented fraction from a logistic regression model.
#' @details Estimates confidence intervals using the delta method on
#' parameters from a generalized linear model with logit link. 
##		Latex code not good here
## \cr \cr  \eqn{RR={{{{\hat{\mu}}}_{2}}}/{{{{\hat{\mu}}}_{1}}}}, where
## \eqn{{\hat{\mu}}_{i}} are 
#' \cr \cr  \eqn{RR={{{\mu}}_{2}}/{{{\mu}}_{1}}}{RR=muhat_2/muhat_1}, where
#' \eqn{{\mu}_{i}}{muhat_i} are 
#' the estimated probabilities from the model.
#' @param fit A \code{\link{glm}} object.
#' @param beta.hat Parameters estimates from a logistic regression with no intercept.
#' @param var.beta.hat Variance-covariance matrix from a logistic regression with no intercept.
#' @param degf Degrees of freedom.
#' @param which Numeric vector indicating which parameters to compare, so that \emph{RR} = compare[2]/compare[1]
#' @param pf Estimate \emph{RR} or its complement \emph{PF}? 
#' @param norm Estimate confidence interval using quantiles of Guassian rather than t distribution quantiles? 
#' @param alpha Complement of the confidence level.
#' @param rnd Number of digits for rounding. Affects display only, not estimates.
#' @return A \code{\link{rror}} object with the following fields.
#'  \item{estimate}{vector with point and interval estimate}
#'  \item{estimator}{either \emph{PF} or \emph{RR}}
#'	\item{mu}{matrix with rows giving probability estimates for each of the groups}
#'  \item{rnd}{how many digits to round the display}
#'  \item{alpha}{complement of confidence level}
#'  \item{norm}{logical indicating Gaussian or t interval}
#'  \item{degf}{degrees of freedom}
#' @export
#' @references No references yet.
#' @author David Siev \email{david.siev@@aphis.usda.gov}
#' @note Level tested: Moderate. (Extensive in S+, not as much in R. They have different \code{glm()} functions.) \cr \cr
#' Call to this function may be one of two formats: (1) specify \code{fit}  or (2) \code{beta.hat}, \code{var.beta.hat}, \code{degf}  \cr \cr
#' \code{RRor(fit, degf = NULL, pf = TRUE, alpha = 0.05, which = c(1, 2), norm = TRUE, rnd = 3)} \cr \cr
#' \code{RRor(beta.hat, var.beta.hat, degf,  pf = TRUE, alpha = 0.05, which = c(1, 2), norm = TRUE, rnd = 3)}
#' @seealso \code{\link{rror}, \link{phiWt}, \link{tauWt}}. See the package vignette for more examples.
#' 
#' @examples
#' bird.fit <- glm(cbind(y, n - y) ~ tx - 1, binomial, bird)
#' RRor(tauWt(bird.fit))
#' 
#' # 95% t intervals on 4 df
#' # 
#' # PF 
#' #     PF     LL     UL 
#' #  0.500 -0.583  0.842 
#' # 
#' #       mu.hat    LL     UL
#' # txcon  0.733 0.943 0.3121
#' # txvac  0.367 0.752 0.0997
#' # 
#' # See the package vignette for more examples

##-----------------------------------------------
## RRor
##-----------------------------------------------

RRor <- function(fit=NULL, beta.hat=NULL, var.beta.hat=NULL, degf=NULL, which = c(1, 2), pf=T, norm = F, alpha = 0.05, rnd=3){
	if(!is.null(fit)){
		beta.hat <- coef(fit)
		var.beta.hat <- summary(fit)$cov.sc
		if(is.null(degf)) degf <- summary(fit)$df.resid
	}
	q <- c(0.5, alpha/2, 1 - alpha/2)
	B <- beta.hat[which]
	b1 <- B[1]
	b2 <- B[2]
	m1 <- 1/(1 + exp( - b1))
	m2 <- 1/(1 + exp( - b2))
	log.r <- log(1 + exp( - b1)) - log(1 + exp( - b2))
	grad.log.r <- c( - exp( - b1) * m1, exp( - b2) * m2)
	var.b <- var.beta.hat[which, which]
	var.log.r <- t(grad.log.r) %*% var.b %*% grad.log.r
	if(norm) {
		int <- exp(log.r + qnorm(q) * sqrt(var.log.r))
		mu <- 1/(1 + exp( - B + matrix(qnorm(q), 2, 3, byrow = T) * sqrt(diag(var.b))))
	}
	else {
		int <- exp(log.r + qt(q, degf) * sqrt(var.log.r))
		mu <- 1/(1 + exp( - B + matrix(qt(q, degf), 2, 3, byrow = T) * sqrt(diag(var.b))))
	}
	dimnames(mu) <- list(names(coef(fit)), c('mu.hat','LL','UL'))
	if(!pf) 
        names(int) <- c("RR", "LL", "UL")
    else{
        int <- 1 - int[c(1,3,2)]
        names(int) <- c("PF", "LL", "UL")
         }
	return(rror$new(estimate = int, estimator = ifelse(pf, 'PF', 'RR'), mu = mu, rnd = rnd, alpha = alpha, 
		norm = norm, degf = degf))
    # out <- list(estimate = int, estimator = ifelse(pf, 'PF', 'RR'), mu = mu, rnd = rnd, alpha = alpha, norm = norm, degf = degf)
    # class(out) <- 'rror'
	# return(out)
}
