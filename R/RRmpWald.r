#' @title Wald confidence intervals for RR from matched pairs
#' @description Estimates confidence intervals for the risk ratio or prevented fraction from matched pairs.
#' @details Estimates confidence intervals for the risk ratio or prevented fraction from matched pairs.
#' The response is the tetranomial vector [11, 12, 21, 22], where the first index is the row and the
#' the second index is the column when displayed as a 2x2 table. Wald type confidence intervals are found by
#' applying the delta method to the multinomial variance. This method fails when there are no responders in one of the treatment groups.
#' \cr \cr Alternative forms of data entry are illustrated by the output, say 
#' \code{Y}, where \code{c(Y$xtable) = Y$freqvec = Y$multvec$Freq}.
#' \cr \cr If RR = 0 (PF = 1), the function will return degenerate interval.
#' @name RRmpWald
#' @param formula Formula of the form \code{y ~ x + cluster(w)}, where y is the 
#' indicator for an individual's positive response, x is a factor with two levels 
#' of treatment, and w identifies the pairs.
#' @param data \code{data.frame} containing variables in formula
#' @param compare Text vector stating the factor levels: compare[1] is the control 
#' or reference group to which compare[2] is compared
#' @param affected Indicator for positive response
#' @param x Alternative data input. Instead of formula and data frame, data may 
#' be input as vector or 2x2 table
#' @param alpha Complement of the confidence level
#' @param pf Estimate \emph{RR} or its complement \emph{PF}?
#' @param tdist Use t distribution? 
#' @param df Degrees of freedom. When NULL, the function will default to 
#' \code{df = N - 2}, where N is the total number of pairs. 
#' @param rnd Number of digits for rounding. Affects display only, not estimates.
#' @return A \code{\link{rrmp}} object with the following fields:
#'  \item{estimate}{vector of point and interval estimates - see details}
#'  \item{estimator}{either \code{"PF"} or \code{"RR"}}
#'  \item{compare}{text vector, same as input}
#'  \item{alpha}{complement of confidence level}
#'  \item{rnd}{how many digits to round the display}
#'  \item{xtable}{data arrayed in 2x2 matrix}
#'  \item{freqvec}{data arrayed a 4-vector}
#'  \item{multvec}{data frame showing the multinomial representation of the data}
#' @export
#' @author David Siev \email{david.siev@@aphis.usda.gov}
#' @note Level tested: Low. \cr\cr Experimental functions for estimating profile likelihood intervals are in the CVBmisc package. \cr \cr
#' Call to this function may be one of two formats: (1) specify \code{data} and \code{formula} or (2) as a vector or 2 x 2 table \code{x} \cr \cr
#' \code{RRmpWald(formula, data, compare = c('con', 'vac'), affected = 1, alpha = 0.05, pf = TRUE, tdist = TRUE, df = NULL, rnd = 3)} \cr \cr
#' \code{RRmpWald(x, compare = c('con', 'vac'), affected = 1, alpha = 0,05, pf = TRUE, tdist = TRUE, df = NULL, rnd = 3)}
#' 
#' @examples
#' RRmpWald(pos ~ tx + cluster(cage), New, compare = c('con', 'vac'))
#'
#' # PF 
#' # 95% interval estimates
#' #
#' #   PF    LL    UL 
#' # 0.550 0.183 0.752 
#' 
#' 
#' RRmpWald(x = c(2, 9, 1, 6))
#' 
#' # PF 
#' # 95% interval estimates
#' 
#' #   PF    LL    UL 
#' # 0.727 0.124 0.915
#' 
RRmpWald <- function(formula = NULL, data = NULL, compare = c('con', 'vac'), 
	affected = 1, x, alpha = 0.05, pf = TRUE, tdist = TRUE, df = NULL, rnd = 3){
	# CI for RR with matched pairs, based on asymptotic normality of log(RR)
	#  and multinomial variance
	# Data entry:
    # formula of the form response ~ treatment + cluster(clustername)
	#   then it will convert data to matrix and vector
	#  if entered as vector (x=) it must be ordered by vac/con pairs: c(11,01,10,00)
	#  if entered as matrix (x=) it must be must be matrix(c(11,01,10,00),2,2,byrow=T)
	multvec <- NULL
	if(!is.null(formula) & !is.null(data)){
        Xx <- .twoby(formula = formula, data = data, compare = compare, affected = affected)
        xtable <- Xx$xtable
		x <- Xx$freqvec
		multvec <- Xx$multvec
	} else if(is.matrix(x)){
		if(!all(dim(x) == 2)){
			stop('Table dimensions must be 2 x 2\n')
		} else {
			xtable <- x
			x <- c(x)
		}
	} else if(is.vector(x)){
		if(length(x)!=4){
			stop('Vector length must be 4\n')
		}
		xtable <- matrix(x,2,2,byrow=T)
	}

	N <- sum(x)
	p <- x/N
	V <- (diag(p) - t(t(p)) %*% t(p)) / N
	p1 <- p[1] + p[2]
	p2 <- p[1] + p[3]
	R <- p2 / p1
	gradR <- c((p[2] - p[3]) / p1^2, -p2 / p1^2, 1 / p1, 0)
	logR <- log(p2) - log(p1)
	gradlogR <- c(1 / p2 - 1 / p1, -1 / p1, 1 / p2, 0)
	varR <- t(gradR) %*% V %*% t(t(gradR))
	varlogR <- t(gradlogR) %*% V %*% t(t(gradlogR))

	if(tdist){
		if(is.null(df))	df <- N - 2
	}
	if(!is.null(df)){
		q <- qt(c(0.5, alpha / 2, 1 - alpha / 2), df)
		what <- paste(100 * (1 - alpha), "% t intervals on ", df, " df\n", sep = "")
	} else {
		q <- qnorm(c(0.5, alpha / 2, 1 - alpha / 2))
		what <- paste(100 * (1 - alpha), "% gaussian interval\n", sep = "")
	}
	
	ci.dl <- exp(logR + q * sqrt(varlogR))
	# use the one based on log(R), unless R=0
	if(R==0){
		ci.d <- R + q * sqrt(varR) 
		ci.dl <- ci.d 
	}

	int <- ci.dl
    if(!pf) {
        names(int) <- c("RR", "LL", "UL")
    } else {
        int <- 1 - int[c(1, 3, 2)]
        names(int) <- c("PF", "LL", "UL")
    }
    
    # out <- list(estimate = int,  estimator = ifelse(pf, 'PF', 'RR'), compare = 
		# compare, alpha = alpha, rnd = rnd, xtable = xtable, freqvec = x, 
		# multvec = multvec)
   # class(out) <- 'rr1'
    # return(out)
	if(is.null(multvec)){

		return(rrmp$new(estimate = int, estimator = ifelse(pf, 'PF', 'RR'), compare = 
			compare, alpha = alpha, rnd = rnd, xtable = xtable, freqvec = x))
	} else {

		return(rrmp$new(estimate = int, estimator = ifelse(pf, 'PF', 'RR'), compare = 
			compare, alpha = alpha, rnd = rnd, xtable = xtable, freqvec = x, multvec =
			multvec))
	}
}
.twoby <- function(formula, data, compare, affected){
	cluster <- function(x) {return(x)}
    this.call <- match.call()
		drop.levels <- function (x){
		for (j in 1:ncol(x)) if (is.factor(x[, j])) 
			x[, j] <- factor(as.character(x[, j]))
		return(x)
		}
	data <- drop.levels(data)
    Terms <- terms(formula,specials = 'cluster', data = data)
	environment(Terms) <- environment()
    A <- model.frame(formula = Terms, data = data)
    dat <- A[, 1]
    group <- A[, 2]
    clusters <- A[, 3]
    tbl <- table(clusters, dat, group)[, 2, 1:2]
	for(i in 1:2){
		tbl[,i] <- ifelse(tbl[,i]==affected,'af','un')
		}
	tbl <- data.frame(tbl)	
	# both levels must be present in both groups
	for(i in 1:2) 
		levels(tbl[,i]) <- c(levels(tbl[,i]),c('af','un')[!c('af','un') %in% levels(tbl[,i])]) 
	xtable <- table(tbl[,compare[2]],tbl[,compare[1]])
	# order table
	xtable <- xtable[c('af','un'),c('af','un')] 
	names(dimnames(xtable)) <- rev(compare)
	dimnames(xtable) <- lapply(dimnames(xtable), function(x){ifelse(x == 'af', 
		'pos', 'neg')})
	freqvec <- c(xtable)
	names(freqvec) <- paste(rep(dimnames(xtable)[[1]], 2), rep(dimnames(xtable)[[2]], 
		c(2, 2)))
	multvec <- as.data.frame(xtable)
	return(list(xtable = xtable, freqvec = freqvec, multvec = multvec))
}

