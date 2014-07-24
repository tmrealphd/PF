#' @title Data class pf
# @name pf-class
#' @aliases pf
#' @rdname pf
#' @section Fields:
#' \itemize{
#' \item[\code{estimator}]{  either \code{"PF"} or \code{"IDR"}}
#' \item[\code{rnd}]{  how many digits to round display}
#' \item[\code{alpha}]{  complement of c.i.}
#' }
#' @seealso \code{\link{rr1}, \link{rrsi}, \link{rrsc}, \link{rrstr}}
#' @export
#' @author Marie Vendettuoli \email{marie.c.vendettuoli@@aphis.usda.gov}
pf <- setRefClass('pf', fields = list(estimator = 'character', rnd = 'numeric',
	alpha = 'numeric'))

	
#' @title Data class rr1
# @name rr1-class
#' @aliases rr1
#' @name rr1
#' @rdname rr1
#' @section Fields:
#' \itemize{
#' \item[\code{estimate}]{  vector with point and interval estimate}
#' \item[\code{estimator}]{  either \code{"PF"} or \code{"IDR"}}
#' \item[\code{y}]{  data vector}
#' \item[\code{rnd}]{  how many digits to round display}
#' \item[\code{alpha}]{  complement of c.i.}
#' }
#' @seealso \code{\link{IDRsc}, \link{RRotsst}, \link{RRtosst}}
#' @export
#' @author Marie Vendettuoli \email{marie.c.vendettuoli@@aphis.usda.gov}
rr1 <- setRefClass('rr1', contains = 'pf', fields = list(estimate = 'numeric',
	y = 'matrix'))

#' @title Data class rror
# @name rror-class
#' @aliases rror
#' @rdname rrorclass
#' @section Fields:
#' \itemize{
#' \item[\code{estimate}]{  vector with point and interval estimate}
#' \item[\code{estimator}]{  either \code{"PF"} or \code{"IDR"}}
#' \item[\code{y}]{  data vector}
#' \item[\code{rnd}]{  how many digits to round display}
#' \item[\code{alpha}]{  complement of c.i.}
#' \item[\code{norm}]{  logical indicating Gaussian or t interval}
#' \item[\code{degf}]{  degrees of freedom}
#' \item[\code{mu}]{  matrix with rows giving probability estimates for each of the groups}
#' }
#' @seealso \code{\link{RRor}}
#' @export
#' @author Marie Vendettuoli \email{marie.c.vendettuoli@@aphis.usda.gov}
rror <- setRefClass('rror', contains = 'rr1', fields = list(norm = 'logical',
	degf = 'numeric', mu = 'matrix'))

#' @title Data class rrsi
# @name rrsi-class
#' @aliases rrsi
#' @rdname rrsi
#' @section Fields:
#' \itemize{
#' \item[\code{y}]{  numeric data vector }
#' \item[\code{k}]{  likelihood ratio criterion}
#' \item[\code{rnd}]{  digits to round display}
#' \item[\code{alpha}]{  complement of c.i.}
#' \item[\code{estimate}]{  vector with point and interval estimate}
#' \item[\code{estimator}]{  either \code{"PF"} or \code{"IDR"}}
#' }
#' @seealso \code{\link{IDRlsi}, \link{RRlsi}}
#' @export
#' @author Marie Vendettuoli \email{marie.c.vendettuoli@@aphis.usda.gov}
rrsi <- setRefClass('rrsi', contains = 'pf', fields = list(y = 'numeric', k = 
	'numeric', estimate = 'numeric'))

#' @title Data class rrmp
# @name rrmp-class
#' @aliases rrmp
#' @rdname rrmp
#' @section Fields:
#' \itemize{
#' \item[\code{estimate}]{  vector with point and interval estimate}
#' \item[\code{estimator}]{  either \code{"PF"} or \code{"IDR"}}
#' \item[\code{y}]{  data vector}
#' \item[\code{rnd}]{  how many digits to round display}
#' \item[\code{alpha}]{  complement of c.i.}
#' \item[\code{xtable}]{  2 x 2 data matrix} 
#' \item[\code{compare}]{  text vector, same as input}
#' \item[\code{frecvec}]{  data arrayed a vector of length 4}
#' \item[\code{multvec}]{  data.frame showing the multinomial representation of the data}
#' }
#' @seealso \code{\link{RRmpWald}}
#' @export
#' @author Marie Vendettuoli \email{marie.c.vendettuoli@@aphis.usda.gov}
rrmp <- setRefClass('rrmp', contains = 'rr1', fields = list(compare = 'character',
	xtable = 'matrix', freqvec = 'numeric', multvec = 'data.frame'))

#' @title Data class rrsc
#' @alias rrsc
#' @rdname rrscclass
#' @section Fields:
#' \itemize{
#' \item[\code{estimate}]{  vector with point and interval estimate}
#' \item[\code{rnd}]{  how many digits to round display}
#' \item[\code{alpha}]{  complement of c.i.}
#' \item[\code{estimator}]{  either \code{"PF"} or \code{"RR"}}
#' \item[\code{y}]{  data vector}
#' }
#' @seealso \code{\link{RRsc}}
#' @export
#' @author Marie Vendettuoli \email{marie.c.vendettuoli@@aphis.usda.gov}
rrsc <- setRefClass('rrsc', contains = 'pf', fields = list(estimate = 'matrix',
	y = 'numeric'))
	
#' @title Data class rrstr
# @name rrstr-class
#' @aliases rrstr
#' @rdname rrstrclass
#' @section Fields:
#' \itemize{
#' \item[\code{estimate}]{  vector with point and interval estimate}
#' \item[\code{rnd}]{  how many digits to round display}
#' \item[\code{alpha}]{  complement of c.i.}
#' \item[\code{estimator}]{  either \code{"PF"} or \code{"RR"}}
#' \item[\code{hom}]{  list of homogeneity statistic, p-value, and degrees of freedom}
#' \item[\code{y}]{  matrix of data}
#' \item[\code{compare}]{  groups compared}
#' }
#' @seealso \code{\link{RRstr}}
#' @export
#' @author Marie Vendettuoli \email{marie.c.vendettuoli@@aphis.usda.gov}
rrstr <- setRefClass('rrstr', contains = 'pf', fields = list(estimate = 'matrix',
	hom = 'list', y = 'matrix', compare = 'character'))
	
