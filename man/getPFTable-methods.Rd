\docType{methods}
\name{getPFTable}
\alias{getPFTable}
\alias{getPFTable,rr1-method}
\alias{getPFTable,rror-method}
\alias{getPFTable,rrsc-method}
\alias{getPFTable,rrsi-method}
\alias{getPFTable,rrstr-method}
\title{Write html and latex table from package PF objects}
\usage{
  getPFTable(x, type = "summary", out = "dev", file = "",
    ...)
}
\arguments{
  \item{x}{object of class \code{rr1, rror, rrsc, rrstr,}
  or \code{rrsi}}

  \item{type}{table type.}

  \item{out}{output type. \code{dev} for eye-readable to
  active device. \code{tex} for latex. \code{texdoc} for
  latex with preamble and footer. \code{doc} for word
  document via R2wd. \code{html} for html file}

  \item{file}{name to save file to. if \code{file == ''},
  display on active device.}

  \item{...}{Additional argument list}
}
\description{
  Write html and latex table from package PF objects
}
\author{
  Marie Vendettuoli
  \email{marie.c.vendettuoli@aphis.usda.gov}
}

