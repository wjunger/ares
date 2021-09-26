\name{setup}
\alias{setup}
\title{Set up ares Environment}
\description{
Set up ares environmental by cacching the data, formatting dates, creating holidays, and subsetting on data.}
\usage{
setup(dataset, date.var, weights = NULL, selection = NULL, 
	date.format = "\%d/\%m/\%Y", weekday.ref = "Sun", holidays = TRUE, ...)
}
\arguments{
  \item{dataset}{the data frame for the analysis}
  \item{date.var}{string with the name of the column that holds the date variable. It should be a text variable or a factor}
  \item{weights}{an optional vector of weights to be used during the estimation process. See Details}
  \item{selection}{indicator variable for sub-setting. It can be either a quoted string with the name of a column in \code{dataset} or a vector}
  \item{date.format}{a string indicating the date format for the horizontal axis. Default is \code{dd/mm/yyyy} See \code{\link{strptime}} for options}
  \item{weekday.ref}{a 3-letter string that match the reference weekday in the U.S. English locale. For example, "Mon" for Monday. Default is "Sun"}
  \item{holidays}{logical indicating if holidays should be created. Default is \code{TRUE}. See Details}
  \item{\dots}{further options for \code{\link{gen_holidays}}}
}
\details{
This function is intended to cache the dataset and set up some environment objects in \R so the rest of the library can work properly. Although most of the functions can work stand alone, some of them depend on the object \code{.ares.active.dataset} and \code{.ares.selection}.

The argument \code{weights} is an optional vector for weighting. Usually, it is a vector of ones if none is specified, but when the user has prior information about the weights of each observations, it can be set accordingly. It is particularly useful when the dataset has suffered some sort of imputation technique.

First, \code{.ares.selection} is assigned to the global environment, then a date formatted variable is created from \code{date.var} and a time index is generated. Factors for days of the week, months, quarters and years are defined. Indicator variables for holidays are created if \code{holidays} is set to \code{TRUE}. The holidays are created by \code{\link{gen_holidays}} relying on the date variable. The date variable is defined as \code{doe} that stands for "date of event". Country-specific holidays can be be generated on set up time for some countries. One can add these country-specific holidays by adding a file named \code{"XX.hol"} to the \code{"/R/library/ares/etc"} directory, where "XX" is a 2-letter abbreviation of the country name, for instance \code{"BR.hol"} stands for the default option Brazil. The holidays file is self documented. Or you can simply send us an email with the holidays and the dates and we add it for you.

Finally, \code{.ares.active.dataset} is assigned to the global environment, every instance of this data frame in the search path is removed and it is reattached.

Attention! The \code{.ares.active.dataset} is replaced without warnings. Any changes made to it during a session that was not saved will be lost.

Note! Be warned that some functions were renamed in this version in order to avoid conflict with S3 methods, namely those starting with either \emph{plot.} or 
\emph{print.}. Please, change your scripts accordingly.

}
\value{
This function invisibly returns \code{.ares.active.dataset} as a data frame. 
}
\references{
Schwartz, J., Spix, C., Touloumi, G. et al. (1996) Methodological issues in studies of air pollution and daily counts of deaths or hospital admissions. \emph{J Epidemiol. Community Health} \bold{50 (suppl 1)}, S12--S18.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

}
\seealso{\code{\link{attach}},\code{\link{search}},\code{\link{as.Date}}}
\examples{
data(admrio)
names(admrio)
setup(admrio,"date")
setup(admrio,"date",selection="warmseason")
}
\keyword{data}
\keyword{datagen}
