#' @title RR likelihood support interval.
#' @descriptionEstimates likelihood support interval for the risk ratio or prevented fraction by the likelihood profile.
#' @details Estimates a likelihood support interval for \emph{RR} or \emph{PF} by the profile likelihood method using the DUD algorithm.
#' \cr \cr Likelihood support intervals are usually formed based on the desired likelihood ratio, often 1/8 or 1/32. Under some conditions the log likelihood ratio
#' may follow the chi square distribution. If so, 
#  Latex version not good
# then \eqn{\alpha=1-{{F}_{{{\chi }^{2}}}}\left( 2\log (k),1 \right)}. 
#' then \eqn{\alpha=1-F(2log(k),1)}, where \eqn{F} is a chi-square CDF. 
#' \code{RRlsi()} will make the conversion from \eqn{\alpha}
#' to \emph{k} if \code{use.alpha = TRUE}.
#' \cr \cr The data may also be a matrix. In that case \code{y} would be entered as \code{matrix(c(y1, n1-y1, y2, n2-y2), 2, 2, byrow = TRUE)}.
# @usage RRlsi(y, k=8, use.alpha = F pf = T)
#' @param y Data vector c(y1, n1, y2, n2) where y are the positives, n are the total, and group 1 is compared to group 2.
#' @param k Likelihood ratio criterion.
#' @param alpha Complement of the confidence level (see details).
#' @param use.alpha Base choice of k on its relationship to alpha? 
#' @param pf Estimate \emph{RR} or its complement \emph{PF}? 
#' @param iter.max Maximum number of iterations
#' @param converge Convergence criterion
#' @param rnd Number of digits for rounding. Affects display onlyRR, not estimates.
#' @param start Optional starting value.
#' @param track Verbose tracking of the iterations? 
#' @param full.track Verbose tracking of the iterations? 
#' @return An object of class \code{\link{rrsi}} with the following fields: \cr
#'  \item{estimate}{matrix of point and interval estimates - see details}
#'  \item{estimator}{either \code{"PF"} or \code{"RR"}}
#'  \item{y}{data vector}
#'  \item{rnd}{how many digits to round the display}
#'  \item{k}{likelihood ratio criterion}
#'  \item{alpha}{complement of confidence level}
#' @export
#' @references Royall R. \emph{Statistical Evidence: A Likelihood Paradigm}. Chapman & Hall, Boca Raton, 1997.  Section 7.6
#' \cr Ralston ML, Jennrich RI, 1978. DUD, A Derivative-Free Algorithm for Nonlinear Least Squares. \emph{Technometrics} 20:7-14. 
#' @author David Siev \email{david.siev@@aphis.usda.gov}
#' @note Level tested: High.
#' @examples
#' RRlsi(c(4,24,12,28))
#' 
#' # 1/8 likelihood support interval for PF 
#' 
#' # corresponds to 95.858% confidence
#' #   (under certain assumptions)
#' 
#' # PF 
#' #     PF     LL     UL 
#' # 0.6111 0.0168 0.8859 

##--------------------------------------------------------------------
## RR likelihood support interval by the profile likelihood method
##--------------------------------------------------------------------
##
## see file RR likelihood functions.r for other versions including conditional likelihood
RRlsi <- function(y, alpha = 0.05, k=8, use.alpha = FALSE, pf=TRUE, iter.max = 50, converge = 1e-006, rnd=3, start = NULL, track = FALSE,full.track = FALSE){
# 9/14/07
# alpha to k: k=exp(qchisq(1-alpha,1)/2)
# k to alpha: alpha=1-pchisq(log(k)*2,1)	
	if(use.alpha)
		k <- exp(qchisq(1-alpha,1)/2)
		else alpha <- 1 - pchisq(log(k) * 2, 1)
		
	# data vector
	if(is.matrix(y))
		y <- c(t(cbind(y[,1],apply(y,1,sum))))
	y1 <- y[1]
	n1 <- y[2]
	y2 <- y[3]
	n2 <- y[4]
	p1 <- y1/n1
	p2 <- y2/n2
	r.hat <- p1/p2

	ll <- function(y, R){
		y1 <- y[1]
		n1 <- y[2]
		y2 <- y[3]
		n2 <- y[4]
		a <- R * (n1 + n2)
		b <-  - (R * (y2 + n1) + y1 + n2)
		cc <- y1 + y2
		det <- sqrt(b^2 - 4. * a * cc)
		p2 <- (- b - det)/(2 * a)
		p1 <- p2 * R
			if(p2 < 1e-4)	p2 <- 1e-004
				else if(p2 > 1 - 1e-4) 	p2 <- 1 - (1e-4)
			if(p1 < 1e-4)	p1 <- 1e-004
				else if(p1 > 1 - 1e-4) 	p1 <- 1 - (1e-4)
		kernel <- y1 * logb(p1) + (n1 - y1) * logb(1 - p1) + y2 * logb(p2) + (n2 - y2) * logb(1 - p2)
		return(kernel)
	}
	ll.max <- y1 * log(p1) + (n1 - y1) * log(1 - p1) + y2 * log(p2) + (n2 - y2) * log(1 - p2)
	end <-  log(1/k) + ll.max
	si <- rep(0,2)
	which <- 1:2
	if(is.null(start)) start <- r.hat*(cbind(c(.4,.5),c(1.5,2)))
	if(r.hat==0){
		si[1] <- 0
		which <- 2
		start[,2] <- start[,2] + c(.6,.9)
	}
	else if(r.hat==1){
		si[2] <- 1
		which <- 1
	}
	for(i in which){
	Q <- c(ll(y,start[1,i]),ll(y,start[2,i]))
	Q <- Q[rev(order(abs(Q-end)))]
	tt <- start[rev(order(abs(Q-end))),i]
	Q2 <- Q[2]
	Q1 <- Q[1]
	t2 <- tt[2]
	t1 <- tt[1]
	iter <- 0
	repeat{
		iter <- iter+1
		t3 <- t2 +(t2-t1)*(end-Q2)/(Q2-Q1)
		if(track) cat('\niter',iter,'r.hat',t2,'LR',exp(Q2-ll.max),'\n')
		if(full.track) cat('t321',t3,t2,t1,'321',ll(y,t3),Q2,Q1,'converge',abs(t3-t2)/t2,'\n')
		if(iter > 1) if(abs((t3-t2)/t2) < converge) break
		t1 <- t2
		t2 <- t3
		Q1 <- Q2
		Q2 <- ll(y,t3)
		if(iter > iter.max) break
		}
	si[i] <- t3
	}

	int <- c(r.hat,si)
	if(!pf) 
        names(int) <- c("RR", "LL", "UL")
    else{
        int <- 1 - int[c(1,3,2)]
        names(int) <- c("PF", "LL", "UL")
        }
    return(rrsi$new(estimate = int, estimator = ifelse(pf, 'PF', 'RR'), y = y, rnd = rnd, k=k, alpha = alpha))

}

