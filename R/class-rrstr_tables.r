#' @export
.summarytable.rrstr.tex <- function(x, file = '', append = F,...){
	cat('%\\usepackage{float, hyperref}\n', file = file, append = append)
	print(xtable(t(as.data.frame(x@hom)), digits = x@rnd,  caption ='Test of homogenity across clusters' ), 
			file = file, append = T, include.colnames = F, hline.after = c(0,3))
	print(xtable(x@estimate, caption = paste(x@estimator, ' ', round(100*(1-x@alpha),x@rnd), 
		'\\% interval estimates', sep = ''), digits = x@rnd), table.placement ='H', file = file, append = T)
}

#' @export
.summarytable.rrstr.texdoc <- function(x, file = '',...){
	cat('\\documentclass[12pt]{article}\n \\usepackage{float, hyperref}\n \\begin{document}\n', file = file)
	.summarytable.rrstr.tex(x = x, file = file, append = T)
	cat('\\end{document}\n', file = file, append = T)
}


#' @export
.summarytable.rrstr.html <- function(x, file = '',...){

	cat('<html>\n<head>\n<link rel="stylesheet" type="text/css" href="W:\\Biometrics\\Section\\Software\\style.css">\n</head>\n<body>\n', file = file)
	print(xtable(t(as.data.frame(x@hom)), digits = x@rnd,  caption ='Test of homogenity across clusters' ), 
			file = file,  type = 'html', 
		sanitize.colnames.function = function(x){paste('<div align = "center">', x, '</div>')}, html.table.attributes = 'width = 300, 
		border = 1', append = T, include.colnames = F, hline.after = c(0,3))
	print(xtable(x@estimate, caption = paste(x@estimator, ' ', round(100*(1-x@alpha),x@rnd), 
		'% interval estimates', sep = '')), file = file, type = 'html', 
		sanitize.colnames.function = function(x){paste('<div align = "center">', x, '</div>')}, html.table.attributes = 'width = 300, 
		border = 1',append = T)
	cat(paste('</body>\n</html>\n', sep =''), file = file, append = T)
}


#' @export
.summarytable.rrstr <- function(x, out = 'dev', file = '',...){
	args <- list(...)
	if(is.null(args$autoformat)){
		autoformat <- 1
	} else {
		autoformat <- args$autoformat
	}
	switch(	out,
			dev = print(x),
			tex = .summarytable.rrstr.tex(x = x, file = file),
			texdoc = .summarytable.rrstr.texdoc(x = x, file = file),
			html = .summarytable.rrstr.html(x = x, file = file),
		)
}

#' @export
.getTable.rrstr <- function(x, type = 'summary', out = 'dev', file = '',...){
	.summarytable.rrstr(x = x, type = type, out = out, file = file,...)
}
