#' Package for PF.
#'
#' Includes functions related to prevented fraction. 
#'
#' \tabular{ll}{
#' Package: \tab pf-package\cr
#' Type: \tab Package\cr
#' Version: \tab 9.5.1\cr
#' Date: \tab 2015-03-31\cr
#' License: \tab MIT\cr
#' LazyLoad: \tab yes\cr
#' LazyData: \tab yes\cr
#' }
#'
#' @name PF-package
#' @aliases PF
#' @docType package
#' @author David Siev \email{David.Siev@@aphis.usda.gov}
#' @examples
#' #---------------------------------------------
#' # Checking PF package
#' #---------------------------------------------
#' example(RRsc)
#' example(RRstr)
#' example(RRmh)
#' example(RRor)
#' example(phiWt)
#' example(tauWt)
#' example(rsbWt)
#' example(rsb)
#' example(RRlsi)
#' example(IDRsc)
#' example(IDRlsi)
#' #---------------------------------------------
#' # The next two take a moment to run
#' #---------------------------------------------
#' example(RRtosst)
#' example(RRotsst)
#' #---------------------------------------------
#' # End examples
#' #---------------------------------------------
#' invisible()
NA

#' @name New
#' @title New dataset
#' @alias New-data
#' @docType data
#' @format a data frame with 52 observations of the following 3 variables, no NAs
#' \describe{
#' \item{cage}{cage ID. 1 - 26}
#' \item{tx}{treatment. one of 'con' or 'vac'}
#' \item{pos}{numeric indicator of positive response. 0 = FALSE or 1 = TRUE}
#' }
#' @references We need some references
#' @keywords datasets
NA

#' @name Set1
#' @title Set1 dataset
#' @alias Set1-data
#' @rdname dfSet1
#' @docType data
#' @format a data.frame with 6 observation of the following 4 variables, no NAs
#' \describe{
#' \item{y}{number positive}
#' \item{n}{total number in group \code{tx} x \code{clus}}
#' \item{tx}{treatment 'vac' or 'con'}
#' \item{clus}{cluster ID}
#' }
#' @references We need some references
#' @keywords datasets
NA

#' @name Table6
#' @title Table6 dataset
#' @alias Table6-data
#' @rdname dfTable6
#' @docType data
#' @format a data.frame with 8 observations of the following 4 variables, no NAs
#' \describe{
#' \item{y}{number positive}
#' \item{n}{total number in group \code{tx} x \code{clus}}
#' \item{tx}{treatment 'a' or 'b'}
#' \item{clus}{cluster ID}
#' }
#' @references Table 1 from Gart (1985) 
#' @keywords datasets
NA

#' @name bird
#' @title bird dataset
#' @alias bird-data
#' @docType data
#' @format a data.frame with 6 observations of the following 4 variables, no NAs
#' \describe{
#' \item{y}{number positive}
#' \item{n}{total number in group \code{tx} x \code{all}}
#' \item{tx}{treatment 'vac' or 'con'}
#' \item{all}{all?}
#' }
#' @references we need some references
#' @keywords datasets
NA

#' @name birdm
#' @title birdm dataset
#' @alias birdm-data
#' @docType data
#' @format a data.frame with 6 observations of the following 4 variables, no NAs
#' \describe{
#' \item{y}{number positive}
#' \item{n}{total number in group \code{tx} x \code{all}}
#' \item{tx}{treatment 'vac' or 'con'}
#' \item{all}{all?}
#' }
#' @references we need some references
#' @keywords datasets
NA

#' @name set1
#' @title set1 dataset
#' @alias set1-data
#' @format a 3 x 4 matrix of data in \code{\link{Set1}}
#' @references we need some references!
#' @keywords datasets
NA

#' @name table6
#' @title table6 dataset
#' @alias table6-data
#' @format matrix for of data in \code{\link{Table6}}
#' @keywords datasets
NA

#' @name rat
#' @title rat dataset
#' @alias rat-data
#' @format a data.frame with 32 observations of the following 3 variables, no NAs
#' \describe{
#' \item{y}{number positive}
#' \item{n}{total number}
#' \item{group}{treatment group: 'control' or 'treated'}
#' }
#' @references Weil's rat data (Table 1 of Rao and Scott)
#' @keywords datasets
NA

