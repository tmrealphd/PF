#' @export
.summarytable.rrsc.tex <- function(x, file = '', append = F,...){
	conf <- paste(x@estimator, " ", round(100*(1-x@alpha),x@rnd),'\\% interval estimates',sep='')
	cat('%\\usepackage{float, hyperref}\n', file = file, append = append)
	print(xtable(x@estimate, caption = conf), table.placement ='H', file = file, append = T)
}

#' @export
.summarytable.rrsc.texdoc <- function(x, file = '',...){
	cat('\\documentclass[12pt]{article}\n \\usepackage{float, hyperref}\n \\begin{document}\n', file = file)
	.summarytable.rrsc.tex(x = x, file = file, append = T)
	cat('\\end{document}\n', file = file, append = T)
}

#' @export
.summarytable.rrsc.html <- function(x, file = '',...){
	conf <- paste(x@estimator, " ", round(100*(1-x@alpha),x@rnd),'% interval estimates',sep='')
	cat('<html>\n<head>\n<link rel="stylesheet" type="text/css" href="W:\\Biometrics\\Section\\Software\\style.css">\n</head>\n<body>\n', file = file)
	print(xtable(x@estimate, caption = conf), file = file, type = 'html', 
		sanitize.colnames.function = function(x){paste('<div align = "center">', x, '</div>')}, html.table.attributes = 'width = 300, 
		border = 1',
		append = T)
	cat(paste('</body>\n</html>\n', sep =''), file = file, append = T)
}

#' @export
.summarytable.rrsc <- function(x, out = 'dev', file = '',...){
	args <- list(...)
	if(is.null(args$autoformat)){
		autoformat <- 1
	} else {
		autoformat <- args$autoformat
	}
	switch(	out,
			dev = print(x),
			tex = .summarytable.rrsc.tex(x = x, file = file),
			texdoc = .summarytable.rrsc.texdoc(x = x, file = file),
			html = .summarytable.rrsc.html(x = x, file = file),
			)
}

#' @export
.getTable.rrsc <- function(x, type = 'summary', out = 'dev', file = '',...){
	.summarytable.rrsc(x = x, type = type, out = out, file = file,...)
}
