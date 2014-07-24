#' Write html and latex table from package PF objects
#' 
#' @param x object of class \code{rr1, rror, rrsc, rrstr,} or \code{rrsi}
#' @param type table type.
#' @param out output type. \code{dev} for eye-readable to active device. \code{tex} for latex. \code{texdoc} for latex with preamble and footer. \code{doc} for word document via R2wd. \code{html} for html file
#' @param file name to save file to. if \code{file == ''}, display on active device.
#' @param ... Additional argument list
#' @export
#' @docType methods
#' @rdname getPFTable-methods
#' @author Marie Vendettuoli \email{marie.c.vendettuoli@@aphis.usda.gov}
getPFTable <- function(x, type = 'summary', out = 'dev', file = '',...){
   attributes(x)
}



#' @rdname getPFTable-methods
#' @aliases getPFTable,rr1-method 
#' @docType methods
setMethod('getPFTable', 'rr1',
	function(x, type = 'summary', out = 'dev', file = '',...){ 
		.getTable.rr1(x = x, type = type, out = out, file = file,...)
	})



#' @rdname getPFTable-methods
#' @aliases getPFTable,rror-method 
#' @docType methods
setMethod('getPFTable', 'rror', 
	function(x, type = 'summary', out = 'dev', file = '',...){ 
		.getTable.rror(x = x, type = type, out = out, file = file,...)
	})



#' @rdname getPFTable-methods
#' @aliases getPFTable,rrsc-method 
#' @docType methods
setMethod('getPFTable', 'rrsc', 
	function(x, type = 'summary', out = 'dev', file = '',...){ 
		.getTable.rrsc(x = x, type = type, out = out, file = file,...)
	})

#' @rdname getPFTable-methods
#' @aliases getPFTable,rrstr-method 
#' @docType methods
setMethod('getPFTable', 'rrstr',
	function(x, type = 'summary', out = 'dev', file = '',...){ 
		.getTable.rrstr(x = x, type = type, out = out, file = file,...)
	})



#' @rdname getPFTable-methods
#' @aliases getPFTable,rrsi-method 
#' @docType methods
setMethod('getPFTable', 'rrsi', 
	function(x, type = 'summary', out = 'dev', file = '',...){ 
		.getTable.rrsi(x = x, type = type, out = out, file = file,...)
	})
