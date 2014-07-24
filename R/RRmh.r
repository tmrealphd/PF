#' @title Mantel-Haenszel method, CI for common RR over strata or clusters with sparse data.
#' @description Estimates confidence intervals for the risk ratio or prevented fraction from clustered or stratified data, using a Mantel-Haenszel estimator for sparse data.
#' @details Based on the Mantel-Haenszel (1959) procedure for sparse data developed by Greenland and Robins (1985).  
#' The confidence limits are based on asymptotic normality of the log(risk ratio).  Agresti and Hartzel (2000) 
#' favor this procedure for small, sparse data sets, but they warn that it is less efficient than maximum likelihood for large data sets.  
#' @param formula Formula of the form \code{cbind(y, n) ~ x + cluster(w)}, where \code{y} is the number positive, \code{n} is the group 
#' size, \code{x} is a factor with two levels of treatment, and \code{w} is a factor indicating the clusters.
#' @param data \code{data.frame} containing variables for formula
#' @param compare Text vector stating the factor levels: \code{compare[1]} is the control or reference group to which \code{compare[2]} is compared
#' @param Y Matrix of data, \eqn{K \times 4}{K x 4}. Each row is a stratum or cluster. The columns are \eqn{y2, n2, y1, n1}, where the y's are 
#' the number of positive in each group, and the n is the total in each group. If data entered by formula and dataframe, \code{Y} is generated automatically.
#' @param pf Estimate \emph{RR} or its complement \emph{PF}? 
#' @param alpha Complement of the confidence level.
#' @param rnd Number of digits for rounding. Affects display only, not estimates.
#' @return An object of class \code{\link{rr1}}  with the following fields.
#'  \item{estimate}{vector of point and interval estimates:  point estimate, lower confidence limit, upper confidence limit}
#'  \item{estimator}{either \code{"PF"} or \code{"RR"}}
#'  \item{y}{y matrix of the data}
#'  \item{compare}{groups compared}
#'  \item{rnd}{how many digits to round the display}
#'  \item{alpha}{complement of confidence level}
#' @export
#' @note If either all y1's or all y2's are zero, a division by zero may occur, and a NaN returned for some values.
#' \cr \cr Level tested:  Low.
#' \cr\cr See vignette \emph{Examples for Stratified Designs} for more examples (enter \code{?PF}).
#' \cr
#' Call to this function may be one of two formats: (1) specify \code{data} and \code{formula} or (2) as a matrix \code{Y} \cr \cr
#' \code{RRmh(formula, data, compare = c('b','a'), pf = TRUE, alpha = 0.05, rnd = 3)} \cr \cr
#' \code{RRmh(Y, pf = TRUE, alpha = 0.05, rnd = 3)}
#' @references Mantel N, Haenszel W, 1959.  Statistical aspects of the analysis of data from retrospective studies of disease.  
#' \emph{Journal of the National Cancer Institute} 22:  719-748.
#' \cr \cr Greenland S, Robins JM, 1985.  Estimation of a common effect parameter from sparse follow-up data.  \emph{Biometrics}  41:  55-68.  Errata, 45:  1323-1324.
#' \cr \cr Agresti A, Hartzel J, 2000.  Strategies for comparing treatments on a binary response with multi-centre data.  \emph{Statistics in Medicine}  19:  1115-1139.
#' \cr Lachin JM, 2000.  \emph{Biostatistical Methods:  The Assessment of Relative Risks} (Wiley, New York), Sec. 4.3.1.
#' @author Christopher Tong \email{Christopher.H.Tong@@aphis.usda.gov}
#' @seealso \code{\link{rr1}}
#' @examples
#' 
#' ## Table 1 from Gart (1985)
#' ##  as data frame
#' RRmh(cbind(y,n) ~ tx + cluster(clus), Table6 , pf = FALSE)
#' 
#' #  RR estimates
#' 
#' # RR 
#' # 95% interval estimates
#' # 
#' #   RR   LL   UL 
#' # 2.67 1.37 5.23 
#' # 
#'
#' ## or as matrix
#' RRmh(Y = table6, pf = FALSE)
#' 
##########################################################################################
#
# Mantel-Haenszel estimate of common risk ratio for K 2x2 tables.
# First draft 5 March 2010.  Updated 11 March 2010, 5 Nov 2010.  
# Revised 28 Dec 2010 to bring input/output format consistent with David Siev's RRstr().
# Note that unlike earlier versions, this version does not check the data types and dimensions
# of the inputs.
# 
##########################################################################################

RRmh <- function(formula = NULL, data = NULL, compare = c('b', 'a'), Y, alpha = 0.05, 
	pf = TRUE, rnd = 3)
{

### Internal function
# 1/18/2012 - this function already exists in the package! don't duplicate as an internal function. mcv
# # # matricize <- function(formula, data,compare=compare){
    # # # # Copied verbatim from David Siev's RRstr.R, version 4.0, in PF package
    # # # assign('cluster',function(x) {return(x)},envir=.GlobalEnv)
    # # # A <- model.frame(formula=formula,data=data)
    # # # A <- data.frame(A[,1],A[,2:3]) # for easier subscripting
    # # # A <- A[order(A[,4],A[,3]),]
    # # # y <- A[,1]
    # # # n <- A[,2]
    # # # x <- as.factor(A[,3])
    # # # if(!any(levels(x)==compare[1]) | !any(levels(x)==compare[2]))
        # # # stop('What is being compared?')
    # # # clus <- A[,4]
    # # # Y1 <- A[x==compare[2],1:2]
    # # # Y2 <- A[x==compare[1],1:2]
    # # # Y <- as.matrix(cbind(Y1,Y2))
    # # # dimnames(Y) <- list(levels(clus),c('y1','n1','y2','n2'))
    # # # return(list(A=A,Y=Y))
# # # }
### End of internal function


# convert data to matrix
    if(!is.null(formula) & !is.null(data)){
        Y <- .matricize(formula = formula, data = data, compare = compare)$Y
    }
# save data and empirical Rs
    Y <- cbind(Y, R.obs = (Y[, 1]/Y[, 2])/(Y[, 3]/Y[, 4]))

### Innards of the algorithm
    zcrit <- qnorm(1 - alpha/2) 
    x1vec <- Y[, 1]
    n1vec <- Y[, 2]
    x2vec <- Y[, 3]
    n2vec <- Y[, 4]
    nk <- n1vec + n2vec

    # Point estimate (Greenland & Robins, eq. 4.  Lachin, eq. 4.17)
    numer <- sum(x1vec * n2vec / (nk))
    denom <- sum(x2vec * n1vec / (nk))
    rr.est <- numer / denom

    # variance of log RR (Greenland & Robins, eq. 13)
    numer <- ( (n1vec * n2vec) * (x1vec + x2vec) - (x1vec * x2vec) * nk )/nk^2
    rk <- x1vec * n2vec / nk
    sk <- x2vec * n1vec / nk
    denom <- sum(rk) * sum(sk)
    var.log.rr <- sum(numer) / denom
    
    # Confidence limits
    rr.ci.lo <- exp(log(rr.est) - zcrit * sqrt(var.log.rr))
    rr.ci.hi <- exp(log(rr.est) + zcrit * sqrt(var.log.rr))

    int <- c(rr.est, rr.ci.lo, rr.ci.hi)

    if(!pf) {
        names(int) <- c("RR", "LL", "UL")
    } else {
        int <- 1 - int[c(1, 3, 2)]
        names(int) <- c("PF", "LL", "UL")
    }
	
    return(rr1$new(estimate = int, estimator = ifelse(pf, 'PF', 'RR'), y = Y, 
		rnd = rnd, alpha = alpha))
    # out <- list(estimate = int,  estimator = ifelse(pf, 'PF', 'RR'), Y = Y, compare = compare, alpha = alpha, rnd=rnd)
    # class(out) <- 'rr1'
    # return(out)
}

