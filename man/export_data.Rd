\name{export_data}
\alias{export_data}
\title{Export Data}
\description{
Export datasets into different formats
}
\usage{
export_data(data, file, text.format = "csv", ...)
}
\arguments{
  \item{data}{a matrix or data frame}
  \item{file}{quoted string with the file name}
  \item{text.format}{rules for text files}
  \item{\dots}{further options for \code{\link{write.table}}}
}
\details{
The function will try to select the proper filter based on the file extension. Valid file types are those supported by the package \pkg{foreign}. There are rules implemented for Stata up to version 9 (*.dta) and xBase files up to version 4 (*.dbf). Excel files (*.xls) are no longer supported. R binary files (*.rda) can also be written directly by the function.

Several text formats can be written by setting \code{text.format} or passing options through \code{\dots}. If none of the above proprietary formats extension is part of the file name, the function will save the data frame in ASCII text format disregarding the extension. Available options for \code{text.format} are \code{"csv"} for comma separated values, \code{"tab"} for columns separated by a tab character and \code{"spc"} for space separated columns and, also, its variants for comma as the decimal separator, \code{"csv2"}, \code{"tab2"} and \code{"spc2"}.
}
\value{
This function does not return a value.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

}
\seealso{\code{\link{attach}},\code{\link{search}}}
\examples{
data(admrio)
export_data(amdrio,"amdrio.dta")
}
\keyword{data}
\keyword{datagen}
