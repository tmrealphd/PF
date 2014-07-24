##### Shared functions for use within package




### this function is used in RRstr and RRmh 
# @title Convert data.frame and formula to a matrix
# @description This function converts a subset of columns in \emph{data} as 
# specified by \emph{formula} into a matrix
# @param formula the formula to use. of format cbind(y, n) ~ tx + cluster(clus)
# @param compare what to compare (length == 2)
# @return list of:
# \itemize{
# \item{\code{A: }}{A data.frame containing only the variables of \emph{formula}}
# \item{\code{Y: }}{a matrix where each compare element is the set of columns (y, n)
# and each unique(clus) is a row}
# }
# @seealso \code{\link{RRmh}, \link{RRstr}}
.matricize <- function(formula, data, compare = compare){

	# 1/18/2012 - added error checking for compare argument. mcv
	# goal: avoid the two following scenarios
	#     if length(compare) < 2: cannot check levels(x) later
	#     if length(compare) > 2: only look at first two compare elements - notify user
	if(length(compare) > 2){
		warning('matricize: length(compare) > 2; only first two elements used')
	} else if(length(compare) < 2){
		stop('matricize: argument compare must have at least two elements')
	}
	
	# 1/18/2012 - added error checking for formula argument. mcv
	# goal: ensure that formula is of meaningful format. if formula is incorrect,
	#      the accessors to A will not work. A better solution would be to access
	# 	   via names.
	if(length(all.vars(formula)) != 4){
		stop('matricize: formula argument must be of format cbind(y, n) ~ tx + cluster(clus)')
	} else {
		if(is.na(pmatch('cbind', strsplit(as.character(formula), '~')[[2]]))){
			stop('matricize: left side of formula argument must be of format cbind(y, n)')
		}
		
		if(length(strsplit(strsplit(as.character(formula), '~')[[3]], '+', fixed = T)[[1]]) != 2){
			stop('matricize: right side of formula argument must be of format tx + cluster(clus)')
		}
	}
	
    cluster <- function(x){return(x)}
	environment(cluster) <- parent.env(environment())
print("ok")
	Terms <- terms(formula, data = data)
	environment(Terms) <- environment()
	A <- model.frame(formula = Terms, data = data)

print("end")	
	A <- data.frame(A[, 1], A[, 2:3]) # for easier subscripting
    A <- A[order(A[, 4], A[, 3]), ]
    y <- A[, 1]
    n <- A[, 2]
    x <- as.factor(A[, 3])
    if(!any(levels(x) == compare[1]) | !any(levels(x) == compare[2])){
		stop('matricize: What is being compared?')
	}
    clus <- A[, 4]
    Y1 <- A[x == compare[2], 1:2]
    Y2 <- A[x == compare[1], 1:2]
    Y <- as.matrix(cbind(Y1, Y2))
    dimnames(Y) <- list(levels(clus), c('y1', 'n1', 'y2', 'n2'))
    return(list(A = A, Y = Y))
}