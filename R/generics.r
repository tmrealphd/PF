#' Print values for PF data obhects.
#' 
#' @param x object of class rr1, rror, rrsi, rrmp, rrstr, rrsc
#' @param ... other arguments not used by this method
#' @rdname print.PF
#' @export
#' @method print rr1
print.rr1 <- function(x,...){
	num.dig <- options()$digits
    options(digits = x$rnd)
	conf <- paste(100*(1 - x$alpha), '%', sep = '')
    cat('\n')
	estimate <- x$estimate
	cat(x$estimator,'\n')
	cat(conf,'interval estimates\n\n')
    print(estimate)
    cat('\n')
    options(digits = num.dig)
    }
	
#' @nord
setMethod('show', 'rr1', function(object){print(object)})

#' @rdname print.PF
#' @method print rror
#' @export
print.rror <- function(x,...){
    num.dig <- options()$digits
    options(digits = x$rnd)
	if(x$norm){
		conf <- paste(100 * (1 - x$alpha), "% gaussian intervals\n", sep = "")
	} else {
		conf <- paste(100 * (1 - x$alpha), "% t intervals on ", x$degf, " df\n", sep = "")
	}
    cat('\n')
	cat(conf)
    cat('\n')
	cat(x$estimator,'\n')
    print(x$estimate)
    cat('\n')
	print(x$mu)
    options(digits = num.dig)
    }
	

#' @nord
setMethod('show', 'rror', function(object){print(object)})



#' @rdname print.PF
#' @method print rrsi
#' @export
print.rrsi <- function(x,...){

	num.dig <- options()$digits
	options(digits = x$rnd)
	support <- paste('1/',round(x$k, x$rnd), ' likelihood support interval for ',
		x$estimator, sep = '')
	conf <- paste(round(100 * (1 - x$alpha), x$rnd), '%', sep = '')
	cat('\n')
	cat(support,'\n')
	 cat('\ncorresponds to', conf, 'confidence\n  (under certain assumptions)\n')
	cat('\n')
	cat(x$estimator,'\n')
	print(x$estimate)
	cat('\n')
	options(digits = num.dig)
}

#' @nord
setMethod('show', 'rrsi', function(object){print(object)})

#' @rdname print.PF
#' @method print rrmp
#' @export
print.rrmp <- function(x,...){
	num.dig <- options()$digits
    options(digits = x$rnd)
	conf <- paste(100 * (1 - x$alpha), '%', sep = '')
    cat('\n')
	estimate <- x$estimate
	cat(x$estimator, '\n')
	cat(conf, 'interval estimates\n\n')
    print(estimate)
    cat('\n')
    options(digits = num.dig)
}

#' @nord
setMethod('show', 'rrmp', function(object){print(object)})


#' @rdname print.PF
#' @method print rrsc
#' @export
print.rrsc <- function(x,...){
    num.dig <- options()$digits
    options(digits = x$rnd)
	conf <- paste(100 * (1 - x$alpha), '%', sep = '')
    cat('\n')
	estimate <- x$estimate[3:5, ]
	cat(x$estimator, '\n')
	cat(conf,'interval estimates\n\n')
    print(estimate)
    cat('\n')
    options(digits = num.dig)
    }
	
#' @nord
setMethod('show', 'rrsc', function(object){print(object)})



#' @rdname print.PF
#' @method print rrstr
#' @export
print.rrstr <- function(x,...){
     num.dig <- options()$digits
    options(digits = x$rnd)
    mle <- x$estimate['mle', 1]
    if(mle == 0 | mle == 1){
        cat("\nHomogeneity test not possible because MLE = ", mle, "\n")
    } else {
		cat("\nTest of homogeneity across clusters\n")
        cat("\nstat\t", x$hom$stat, "\ndf\t", x$hom$df, "\np\t", x$hom$p, "\n")
        if(x$hom$p <= x$alpha){
            cat("\n\tHeterogeneity may be present\n")
        }
	}
	conf <- paste(100 * (1 - x$alpha), '%', sep = '')
    cat('\n')
	cat(x$estimator, '\n')
	cat(conf, 'interval estimates\n\n')
    print(x$estimate)
    cat('\n')
    options(digits = num.dig)
}
	

#' @nord
setMethod('show', 'rrstr', function(object){print(object)})
