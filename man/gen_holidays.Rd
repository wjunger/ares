\name{gen_holidays}
\alias{gen_holidays}
\title{Generate Holidays}
\description{
Generate indicator variables for holidays and other special dates
}
\usage{
gen_holidays(date, holidays = NULL, dates = NULL, 
	selection = TRUE, country = NULL)
}
\arguments{
  \item{date}{date formatted vector. See \code{\link{as.Date}}}
  \item{holidays}{vector of strings for the holidays names}
  \item{dates}{vector of strings for the holidays dates formatted as \code{"dd/mm"} or \code{"dd/mm/yyyy"}. See Details}
  \item{selection}{a logical indicating whether \code{.ares.selection} should honoured}
  \item{country}{a string indicating the country rules. Default is \code{NULL}. See Details}
}
\details{
This function scans \code{date} looking for matches to each element of \code{dates} and if it finds any 1 is assigned to the variable 0 otherwise. The column is named according to \code{holidays}. If both \code{holidays} and \code{date} are omitted and \code{country} is set to a valid rule, country-specific holidays are included in the matrix. If both \code{.holidays} and \code{.dates} vectors are set in the global environment they are used instead of country rules (that are overridden) to flag the fixed holidays. Moving holidays based on Easter are computed by the function \code{\link{moving.holidays}}.
 
Automatically created holidays are Christmas, New year, Carnaval, Passion, Easter and Corpus Christi. If both \code{holidays} and \code{date} are supplied by the function caller then the resulting matrix will contain only those holidays and these automatically created holidays will be omitted. Thus, one can use this function for both automatically create holidays and/or some special days indicators. The vector \code{dates} can contain dates formatted as both \code{"dd/mm"} and \code{"dd/mm/yyyy"}. The first format will flag all occurrences of that day and month in \code{date} and the second, obviously, will flag the only one occurrence.

Other country rules can be added to this function in the future. If you have a list of holidays for you country, send us a note.
}
\value{
A matrix of indicator variables for the holidays with as many rows as \code{date}.
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
# some mexican holidays plus moving holidays
holidays <- c("constitucion","cincomayo","independencia")
dates <- c("05/02","05/05","16/09")
doe <- .aresEnv$ares.active.dataset$doe
mx.holidays <- gen_holidays(doe)

# special days
strikes <- gen_holidays(doe,holidays=c("transportation","hospitals"),
	dates=c("02/04/2001","21/05/2001"))
}
\keyword{data}
\keyword{datagen}
