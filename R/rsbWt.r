#' @title Rao-Scott weighting. 
#' @description Rao-Scott weighting of clustered binomial observations.
#' @details Estimates the cluster design effect \eqn{{d}_{i}}{d_i} as the variance inflation due to clustering by the method of Rao and Scott. 
#' Observations are then weighted by the inverse of the \eqn{{d}_{i}}{d_i}. 
#' @param fit A \code{\link{glm}} object.
#' @param subset.factor Factor for estimating the weights by subset. 
#' @param fit.only Return only the new fit? If FALSE, also returns the weights and phi estimates.
#' @return A list with the following elements.
#'  \item{fit}{the new model fit, updated by the estimated weights}
#'  \item{weights}{vector of weights}
#'	\item{d}{vector of \eqn{{d}_{i}}{d_i} estimates}
#' @export
#' @references Rao JNK, Scott AJ, 1992. A simple method for the analysis of clustered binary data. \emph{Biometrics} 48:577-585.
#' @author David Siev \email{david.siev@@aphis.usda.gov}
#' @note Level tested: Moderate. 
#' @seealso \code{\link{RRor}, \link{rsb}}. See the package vignette for more examples.
#' 
#' @examples
#' birdm.fit <- glm(cbind(y,n-y)~tx-1,binomial,birdm)
#' RRor(rsbWt(birdm.fit))
#' # 
#' # 95% t intervals on 4 df
#' #
#' # PF 
#' #     PF     LL     UL 
#' #  0.479 -1.061  0.868 
#' #
#' #       mu.hat    LL     UL
#' # txcon  0.768 0.968 0.2659
#' # txvac  0.400 0.848 0.0737
#' #
#' # See the package vignette for more examples
##-------------------------------
## rsbWt uses rsb to refit model
##-------------------------------
rsbWt <- function(fit=NULL, subset.factor=NULL, fit.only = T){
	fit <- update(fit, x = T, y = T)
	x <- fit$x
	yovern <- fit$y
	n <- fit$prior.weights
	if(is.null(n))
		n <- rep(1, length(fit$y))
	y <- yovern * n
	if(is.null(subset.factor))
		subset.factor <- factor(rep('all',length(y)))
	rsbdw <- rsb(y, n, subset.factor)
	w <- rsbdw$w
	d <- rsbdw$d
	newfit <- update(fit,weights = w)
	if(fit.only) out <- newfit
		else out <- list(fit = newfit, weights = w, d = d)
	return(out)
	}

#' @title Rao-Scott weights.
#' @details Estimates the cluster design effect \eqn{{d}_{i}}{d_i} as the variance inflation due to clustering by the method of Rao and Scott. 
#' \code{rsb} estimates the \eqn{{d}_{i}}{d_i} for use by \code{rsbWt} or other functions.
#' @description Rao-Scott weights. 
#' @param y Number positive.
#' @param n Total number.
#' @param id Factor for estimating the weights by subset.
#' @return A list with the following elements.
#'  \item{w}{vector of weights}
#'	\item{d}{vector of \eqn{{d}_{i}}{d_i} estimates}
#' @export
#' @references Rao JNK, Scott AJ, 1992. A simple method for the analysis of clustered binary data. \emph{Biometrics} 48:577-585.
#' @author David Siev \email{david.siev@@aphis.usda.gov}
#' @note Level tested: Moderate. 
#' @seealso \code{\link{rsbWt}}. See the package vignette for more examples.
#' @examples
#' # Weil's rat data (Table 1 of Rao and Scott)
#' rsb(rat$y, rat$n, rat$group)$d
#' #  control  treated 
#' # 1.232495 3.952861 
#-------------------------------
# rsb returns d's and weights
#-------------------------------
rsb <- function(y, n, id=NULL){
	# Rao-Scott design effect weights
	# Biometrics 48:577-585
	# based on raoscott.bin function
	if(is.null(id)) id <- rep('all',length(y))
	y.i <- tapply(y, id, sum)
	n.i <- tapply(n, id, sum)
	m.i <- c(table(id))
	p.i <- y.i/n.i
	r.ij <- y - n * p.i[tapply(y, id)]
	v.i <- (tapply(r.ij^2, id, sum) * m.i)/((m.i - 1) * n.i^2)
	d.i <- (n.i * v.i)/(p.i * (1 - p.i))
	w <- 1 / d.i[tapply(y, id)] 
	y.tilde <- y * w
	n.tilde <- n * w
	return(list(weights=w, d=d.i))
	}
	
