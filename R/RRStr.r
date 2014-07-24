#' @title Gart-Nam method, CI for common RR over strata or clusters.
#' @description Estimates confidence intervals for the risk ratio or prevented fraction from clustered or stratified data.
#' @details Uses the DUD algorithm to estimate confidence intervals by the method of Gart.
#' @param formula Formula of the form \code{cbind(y, n) ~ x + cluster(w)}, where y is the number positive, n is the group size, x is a factor with two levels of treatment, and w is a factor indicating the clusters.
#' @param data data.frame containing variables of formula
#' @param compare Text vector stating the factor levels: compare[1] is the control or reference group to which compare[2] is compared
#' @param Y Matrix of data. Each row is a stratum or cluster. The columns are y2, n2, y1, n1. If data entered by formula and dataframe, Y is generated automatically.
#' @param pf Estimate \emph{RR} or its complement \emph{PF}? 
#' @param alpha Size of the homogeneity test and complement of the confidence level.
#' @param trace.it verbose tracking of the iterations? 
#' @param iter.max Maximum number of iterations
#' @param converge Convergence criterion
#' @param rnd Number of digits for rounding. Affects display only, not estimates.
#' @param multiplier internal control parameter for algorithm
#' @param divider internal control parameter for algorithm
#' @return A \code{\link{rrstr}} object with the following fields:
#'  \item{estimate}{matrix of point and interval estimates - starting value, MLE, and skewness corrected}
#'  \item{hom}{list of homogeneity statistic, p-value, and degrees of freedom}
#'  \item{estimator}{either \code{"PF"} or \code{"RR"}}
#'  \item{y}{Y matrix of the data}
#'  \item{compare}{groups compared}
#'  \item{rnd}{how many digits to round the display}
#'  \item{alpha}{size of test; complement of confidence level}
#' @export
#' @note See vignette \emph{Examples for Stratified Designs} for more examples (enter \code{?PF}).
#' @references Gart JJ, 1985. Approximate tests and interval estimation of the common relative risk in the combination of \eqn{2 x 2} tables. \emph{Biometrika} 72:673-677.
#' \cr Gart JJ, Nam J, 1988. Approximate interval estimation of the ratio of binomial parameters: a review and corrections for skewness. \emph{Biometrics} 44:323-338.
#' \cr Ralston ML, Jennrich RI, 1978. DUD, A Derivative-Free Algorithm for Nonlinear Least Squares. \emph{Technometrics} 20:7-14. 
#' @author David Siev \email{david.siev@@aphis.usda.gov}
#' @note Level tested: High. \cr \cr
#' Call to this function may be one of two formats: (1) specify \code{data} and \code{formula} or (2) as a matrix \code{Y} \cr \cr
#' \code{RRstr(formula, data, compare = c('b','a'), pf = TRUE, alpha = 0.05, trace.it = FALSE, iter.max = 24, converge = 1e-6, rnd = 3, multiplier = 0.7, divider = 1.1)} \cr \cr
#' \code{RRstr(Y, compare = c('b','a'), pf = TRUE, alpha = 0.05, trace.it = FALSE, iter.max = 24, converge = 1e-6, rnd = 3, multiplier = 0.7, divider = 1.1)}
#' @seealso \code{\link{rrstr}}
#' @examples
#' ## Table 1 from Gart (1985)
#' ##  as data frame
#' RRstr(cbind(y,n) ~ tx + cluster(clus), Table6 , pf = FALSE)
#' 
#' # Test of homogeneity across clusters
#' 
#' # stat     0.954 
#' # df       3 
#' # p        0.812 
#' 
#' #  RR estimates
#' 
#' # 		        RR   LL   UL
#' # starting   2.66 1.37 5.18
#' # mle        2.65 1.39 5.03
#' # skew corr  2.65 1.31 5.08
#'
#' ## or as matrix
#' RRstr(Y = table6, pf = FALSE)
#' 
##--------------------------------------------------------------------
## RRstr function
##--------------------------------------------------------------------
RRstr <- function(formula = NULL, data = NULL, compare = c('b', 'a'), Y, alpha = 0.05, 
	pf = TRUE, trace.it = FALSE, iter.max = 24, converge = 1e-6, rnd = 3, multiplier = 0.7, 
	divider = 1.1){
# stratified risk ratio (Gart and Nam, 1988)
#  coded by D. Siev 1993
#  modified by Siev 3/14/2010 to fix bug
#  modified by Siev 5/12/2010 to stabilize algorithm
#   5/17/2010 moved optimization into separate internal function
#   9/24/2010 added formula-dataframe data entry
#	11/25/2010 minor change in convergence

#---------------------------------------
# define internal functions:
#  u.p, zi.phi, zis.phi, root, rr.opt, matricize
 
    u.p <- function(p1, p2, n1, n2)
    (1 - p1)/(n1 * p1) + (1 - p2)/(n2 * p2)

    zi.phi <- function(phi, Y, u.p, root, za){
    # for score interval in RRstr
        zi <- ui <- 0
        for(i in 1:nrow(Y)) {
            y1 <- Y[i, 1]
            n1 <- Y[i, 2]
            y2 <- Y[i, 3]
            n2 <- Y[i, 4]
            p2 <- root(y1, y2, n1, n2, phi)
            p1 <- p2 * phi
            u <- u.p(p1, p2, n1, n2)
            u[u<0] <- NA
            zz <- ((y1 - n1 * p1)/(1 - p1))
            z <- zz * sqrt(u)
            zd <- abs(z - za)
            zd[is.na(zd)] <- 2 * zd[!is.na(zd)]
            zi <- zi + unique(zz[zd == min(zd)])
            ui <- ui + 1/unique(u[zd == min(zd)])
        }
        return(zi/sqrt(ui))
    }

    zis.phi <- function(phi, Y, u.p, root, za){
    # for skewness-corrected interval in RRstr
        zi <- ui <- gi <- 0
        for(i in 1:nrow(Y)) {
            y1 <- Y[i, 1]
            n1 <- Y[i, 2]
            y2 <- Y[i, 3]
            n2 <- Y[i, 4]
            p2 <- root(y1, y2, n1, n2, phi)
            p1 <- p2 * phi
            q1 <- 1 - p1
            q2 <- 1 - p2
            u <- u.p(p1, p2, n1, n2)
            u[u<0] <- NA
            zz <- ((y1 - n1 * p1)/(1 - p1))
            z.ph <- zz * sqrt(u)
            g <- (q1 * (q1 - p1))/(n1 * p1)^2 - (q2 * (q2 - p2))/(n2 * p2)^2
            g.ph <- g/u^1.5
            z.s <- z.ph - (g.ph * (za^2 - 1))/6
            zd <- abs(z.s - za)
            zd[is.na(zd)] <- 2 * zd[!is.na(zd)]
            zi <- zi + unique(zz[zd == min(zd)])
            ui <- ui + 1/unique(u[zd == min(zd)])
            gi <- gi + unique(g[zd == min(zd)])/unique(u[zd == min(zd)])^3
        }
        return(zi/sqrt(ui) - (gi * (za^2 - 1))/(6 * ui^1.5))
    }

    root <- function(y1, y2, n1, n2, phi){
		# in RRstr
        a <- phi * (n1 + n2)
        b <-  - (phi * (y2 + n1) + y1 + n2)
        cc <- y1 + y2
        det <- sqrt(b^2 - 4 * a * cc)
        r1 <- ( - b + det)/(2 * a)
        r2 <- ( - b - det)/(2 * a)
        r1 <- round(r1, 8) ## round for comparison
        r2 <- round(r2, 8) ##
        if(r1 < 0 || r1 > 1) r1 <- r2
        if(r2 < 0 || r2 > 1) r2 <- r1
        return(c(r1, r2))
    }

    rr.opt <- function(z.phi, phi,za, divider, trace.it, u.p, root){
    # optimizer function for RRstr
        zz <- c(z.phi(phi[1], Y, u.p, root, za), z.phi(phi[2], Y, u.p, root, za))
        if(abs(za - zz[1]) > abs(za - zz[2]))  phi <- rev(phi)
        phi.new <- phi[1]
        phi.old <- phi[2]
        z.old <- z.phi(phi.old, Y, u.p, root, za)
        if(trace.it) cat('\n\nR start', phi, '\n')
        iter <- 0
        repeat {
            iter <- iter + 1
            z.new <- z.phi(phi.new, Y, u.p, root, za)
            if(is.na(z.new)) f <- 1 ## adjust step
            while(is.na(z.new)){
                phi.new <- phi.old + (phi.new - phi.old) / f
                z.new <- z.phi(phi.new, Y, u.p, root, za)
                f <- f * (1/divider) 
                if(trace.it) cat('Step adjustment =', 1/f, '\n')
            }

            phi <- exp(log(phi.old) + log(phi.new/phi.old) * ((za - z.old)/(z.new - z.old)))
            phi.old <- phi.new
            z.old <- z.new
            phi.new <- phi

            if(trace.it)
                cat("iteration", iter, "  z", z.new, "phi", phi.new, "\n")
            if(abs(za - z.new) < converge) break
            if(iter == iter.max){ # no convergence 
                cat('\nIteration limit reached without convergence\n')
                break
            }
        } # end repeat
    return(phi.new)
	}

# 1/18/2012 moved to aab.r this is a shared function do not duplicate as an 
# internal function
# # matricize <- function(formula, data, compare = compare){
    # # assign('cluster',function(x) {return(x)}, envir = .GlobalEnv)
    # # A <- model.frame(formula = formula, data = data)
    # # A <- data.frame(A[,1],A[,2:3]) # for easier subscripting
    # # A <- A[order(A[,4],A[,3]),]
    # # y <- A[,1]
    # # n <- A[,2]
    # # x <- as.factor(A[,3])
    # # if(!any(levels(x)==compare[1]) | !any(levels(x)==compare[2]))
        # # stop('What is being compared?')
    # # clus <- A[,4]
    # # Y1 <- A[x==compare[2],1:2]
    # # Y2 <- A[x==compare[1],1:2]
    # # Y <- as.matrix(cbind(Y1,Y2))
    # # dimnames(Y) <- list(levels(clus),c('y1','n1','y2','n2'))
    # # return(list(A=A,Y=Y))
# # }

#---------------------------------------
# end internal function definitions
#---------------------------------------
    this.call <- match.call()

# convert to matrix

    if(!is.null(formula) & !is.null(data)){
        Y <- .matricize(formula = formula, data = data, compare = compare)$Y
    }

# save data and empirical Rs
    Y <- cbind(Y, R.obs = (Y[, 1]/Y[, 2])/(Y[, 3]/Y[, 4]))

    zv <- qnorm(c(alpha / 2, 1 - alpha / 2))
    numer <- denom <- uu <- 0
    for(i in 1:nrow(Y)) {
        y1 <- Y[i, 1]
        n1 <- Y[i, 2]
        y2 <- Y[i, 3]
        n2 <- Y[i, 4]
        p1 <- y1 / n1
        p2 <- y2 / n2
        numer <- numer + (n2 * y1)/max((n1 + n2 - y1 - y2), 1) # added max
        denom <- denom + (n1 * y2)/max((n1 + n2 - y1 - y2), 1) #
        uu <- uu + ifelse(u.p(p1, p2, n1, n2) > 0, 1 / u.p(p1, p2, n1, n2), 0)
    } # end for i

    phi <- Phi <- numer/denom
    v <- sqrt(uu)
    int <- sort(exp(log(phi) + (log(phi) * zv)/v))
    int[int < 0] <- 0

    # if all zeros or ones estimate only one end of interval
    which <- 1:2 
    score <- rep(0, 2)
    if(numer == 0){
        int[1] <- score[1] <- 0
        int[2] <- 1
        which <- 2 
    }
    if(denom == 0){
        int[2] <- score[2] <- 1
        int[1] <- 0
        which <- 1 
    }
    if(trace.it) {
        cat("\nstarting estimates:")
        print(int)
    }

#---------------------------------------
# get MLE of R
# and test for heterogeneity
    if(Phi == 0 | Phi == 1){
        Phi.ML <- Phi
        hom <- paste('MLE = ',Phi.ML,', Homogeneity test not possible', sep = '')
    } else{
        za <-  0
        phi <- c(Phi, multiplier * Phi)
        if(trace.it) cat('\nMLE')
        phi.new <- rr.opt(z.phi = zi.phi, phi = phi, za = za, divider = divider,
			trace.it = trace.it, u.p = u.p, root = root)
        Phi.ML <- phi.new
            # test for homogeneity of phi's
        test <- rep(0, nrow(Y))
        for(i in 1:nrow(Y)){
            test[i] <- zi.phi(Phi.ML, t(as.matrix(Y[i, ])), u.p, root, za)^2
		}
        hom.test <- sum(test)
        df.test <- nrow(Y) - 1
        p.test <- 1 - pchisq(hom.test, df.test)
        hom <- list(stat = hom.test, df = df.test, p = p.test)
    } 
#---------------------------------------

# Score intervals
    for(k in which) {
        if(trace.it) cat('\nScore',switch(k, 'lower', 'upper'))
        za <-  -zv[k]
        phi <- c(int[k], multiplier * int[k])
        phi.new <- rr.opt(z.phi = zi.phi, phi = phi, za = za, divider = divider,
			trace.it = trace.it, u.p = u.p, root = root)
        score[k] <- phi.new
    }
    int <- rbind(int, score)

# Skewness-corrected intervals
    for(k in which) {
        if(trace.it) cat('\nSkewness-corrected',switch(k,'lower','upper'))
        za <-  -zv[k]
        phi <- c(int[1, k], multiplier * int[1, k])
        phi.new <- rr.opt(z.phi = zis.phi, phi = phi, za = za, divider = divider,
			trace.it = trace.it, u.p = u.p, root = root)
        score[k] <- phi.new
    }   # end for k

    int <- rbind(int, score)
    int <- cbind(c(Phi, rep(Phi.ML, nrow(int) - 1)), int)
    if(trace.it) cat("\n\n")
    if(!pf){
		dimnames(int) <- list(c("starting", "mle", "skew corr"), c("RR", "LL", "UL"))
	} else {
        int <- 1 - int[, c(1, 3, 2)]
        dimnames(int) <- list(c("starting", "mle", "skew corr"), c("PF", "LL", "UL"))
    }
    return(rrstr$new(estimate = int, hom = hom, estimator = ifelse(pf, 'PF', 'RR'), 
		y = Y, compare = compare, rnd = rnd, alpha = alpha))
    # out <- list(estimate = int, hom = hom, estimator = ifelse(pf, 'PF', 'RR'), Y = Y, compare = compare, rnd = rnd, alpha = alpha)
    # class(out) <- 'rrstr'
    # return(out)
}

