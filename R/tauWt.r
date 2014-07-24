#' @title Binomial dispersion: intra-cluster correlation parameter. 
#' @description MME estimates of binomial dispersion parameter tau (intra-cluster correlation).
#' @details Estimates binomial dispersion parameter \eqn{\tau} by the method of moments. Iteratively refits the model by the Williams 
#' procedure, weighting the observations by \eqn{1/\phi_{ij}}{1/\phi_ij},
#' where \eqn{\phi_{ij}=1+\tau _j(n_{ij}-1)}{\phi_ij=1+\tau_j(n_ij - 1)}, 
#' \eqn{j} indexes the subsets, and \eqn{i} indexes the observations. 
#' @param fit A \code{\link{glm}} object.
#' @param subset.factor Factor for estimating tau by subset. 
#' @param fit.only Return only the final fit?  If FALSE, also returns the weights and tau estimates.
#' @param iter.max Maximum number of iterations.
#' @param converge Convergence criterion: difference between model degrees of freedom and Pearson's chi-square. Default 1e-6.
#' @param trace.it Display print statments indicating progress
#' @return A list with the following elements.
#'  \item{fit}{the new model fit, updated by the estimated weights}
#'  \item{weights}{vector of weights}
#'	\item{phi}{vector of phi estimates}
#' @export
#' @references Williams DA, 1982. Extra-binomial variation in logistic linear models. \emph{Applied Statistics} 31:144-148.
#' \cr Wedderburn RWM, 1974. Quasi-likelihood functions, generalized linear models, and the Gauss-Newton method. \emph{Biometrika} 61:439-447.
#' @author David Siev \email{david.siev@@aphis.usda.gov}
#' @note Level tested: Moderate. (Extensive in S+, not as much in R. They have different \code{glm()} functions.)
#' @seealso \code{\link{phiWt}}, \code{\link{RRor}}, see package vignette for more examples
#' 
#' @examples
#' birdm.fit <- glm(cbind(y,n-y)~tx-1, binomial, birdm)
#' RRor(tauWt(birdm.fit))
#'  
#' # 95% t intervals on 4 df
#' # 
#' # PF 
#' #     PF     LL     UL 
#' #  0.489 -0.578  0.835 
#' # 
#' #       mu.hat    LL    UL
#' # txcon  0.737 0.944 0.320
#' # txvac  0.376 0.758 0.104
#' #
#' # See the package vignette for more examples

# binomial family only
# any link
tauWt <- function(fit,subset.factor=NULL, fit.only = TRUE, iter.max = 12, converge = 1e-6, trace.it = FALSE){
	# define functions Tau and tauOptim
	Tau <- function(fit, x, n, tau.hat, iter){
		pcs <- sum(resid(fit, type = "pearson")^2)
		degf <- summary(fit)$df.resid
		mu.hat <- fit$fitted
		V.beta <- summary(fit)$cov.sc
		V.eta <- x %*% V.beta %*% t(x)
		d <- diag(V.eta)
		w <- 1 / (1 + (n-1) * tau.hat)
		v <- n * mu.hat * (1 - mu.hat)
		wvd <- w * (1 - w * v * d)
		denominator <- sum((n-1) * wvd)
		if(iter==1)
			numerator <- pcs - degf
			else
			numerator <- pcs - sum(wvd)
		tau.new <- numerator / denominator
		return(tau.new)
		}
		
	tauOptim <- function(fit, x, n){
		tau.hat <- 0
		w <- rep(1, length(n))
		pcs <- sum(resid(fit, type = "pearson")^2)
		df <- fit$df.residual
		iter <- 0
		if(trace.it) cat("\nCycle", iter, "   tau.hat =", tau.hat, "PCS =", pcs, "deviance =", fit$deviance)
		repeat {
			iter <- iter + 1
			if(iter > iter.max)
				break
			tau.hat <- Tau(fit, x, n, tau.hat, iter)
			w <- 1/(1 + tau.hat * (n - 1))
			fit <- update(fit, weights = w, maxit = iter.max/2)
			pcs <- sum(resid(fit, type = "pearson")^2)
			if(trace.it) cat("\nCycle", iter, "   tau.hat =", tau.hat, "PCS =", pcs, "deviance =", round(fit$deviance, 3))
			if(pcs < (df + converge))
				if(iter > 1)
					break
			}
		if(trace.it) cat("\ntau =", tau.hat, "\n")
		return(list(w=w,tau.hat=tau.hat))
		}
		
	fit <- update(fit, x = T, y = T)
	x <- fit$x
	yovern <- fit$y
	n <- fit$prior.weights
	if(is.null(n))
		n <- rep(1, length(fit$y))
	y <- yovern * n
	nminusy <- n - y

	if(is.null(subset.factor)){
		wTau <- tauOptim(fit, x, n)
		w <- wTau$w
		tau.hat <- wTau$tau.hat
		}
	else{
		w <- rep(NA, length(n))
		tau.hat <- rep(NA, length(levels(subset.factor)))
		names(tau.hat) <- levels(subset.factor)
		flink <- fit$family$link ##
		for(lev in levels(subset.factor)){
			if(trace.it) cat('\n',lev,sep='')
			subdat <- data.frame(
				yi = y[subset.factor==lev],
				nminusyi = nminusy[subset.factor==lev]
				)
			subfit <- glm(cbind(yi,nminusyi)~1,binomial(link=flink),data=subdat)
			nobs <- sum(subset.factor==lev)
			xi <- matrix(1,nobs,1)
			wTau <- tauOptim(subfit, xi, n[subset.factor==lev])
			w[subset.factor==lev] <- wTau$w
			tau.hat[lev] <- wTau$tau.hat
			}
		}
	newfit <- update(fit,weights = w)
	if(fit.only) out <- newfit
		else out <- list(fit = newfit, weights = w, tau = tau.hat)
	return(out)
}
