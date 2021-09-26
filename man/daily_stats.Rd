\name{daily_stats}
\alias{daily_stats}
\title{Daily Statistics}
\description{
Compute daily basis statistics from hourly measures
}
\usage{
daily_stats(dataset, parameter, first.column = 2, date = TRUE,
	samples = 24, statistic = "mean", daylight=c(6,19), date.format = "\%d/\%m/\%Y")
}
\arguments{
  \item{dataset}{input dataset}
  \item{parameter}{quoted string indicating the environmental parameter rules}
  \item{first.column}{index of the first column of actual data. The first column is usually allocated for dates, so default is 2}
  \item{date}{a logical indicating whether the resulting dataset should contain a date column}
  \item{samples}{an integer indicating the number of samples per day}
  \item{statistic}{a quoted string indicating the aggregation function. Default is \code{"mean"}. See Details}
  \item{daylight}{a vector indicating the daylight hours. Default is \code{c(,6,19)}}
  \item{date.format}{a string indicating the date format for the horizontal axis. Default is \code{dd/mm/yyyy}. See \code{\link{strptime}} for options}
}
\details{
This is an environmental-specific function. Suppose, for instance, that during the day 24 measurements of a pollutant concentration are taken and one wants to compute daily statistics about the monitored parameter. Besides, when \code{parameter} is set to some \dQuote{known} air pollutant, the aggregation is computed according to specific predefined rules.

The available statistics are: \code{"mean"}, \code{"min"}, and \code{"max"}. Whilst \emph{min} and \emph{max} are the trivial maximum and minimum respectively, \emph{mean} are computed using different rules for each pollutant.

The implemented rules for the mean of some pollutants are the following:

\code{TEMP} or \code{HUMID}: ordinary daily mean are computed;

\code{PM10} or \code{SO2}: daily mean are computed if at least three quarters of hourly observations is not missing;

\code{NO2} or \code{O3}: daily mean are computed if at least three quarters of hourly observations is not missing and at least one quarter lies between 6 a.m and 6 p.m.;

\code{CO}: if parameter is set to \code{CO}, statistic is ignored and the maximum 8-hour running mean is computed.

Attention! Be sure that all columns on the dataset share the same rule for statistic computation, otherwise it will produce wrong values.

Both \code{parameter} and \code{statistic} are not case sensitive.
}
\value{
A data frame with the same columns there are in \code{dataset} plus the date column.
}
\references{
Schwartz, J., Spix, C., Touloumi, G. et al. (1996) Methodological issues in studies of air pollution and daily counts of deaths or hospital admissions. \emph{J Epidemiol. Community Health} \bold{50 (suppl 1)}, S12--S18.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

}
\seealso{\code{\link{import_data}},\code{\link{export_data}}}
\examples{
## Not run:
# pm.mean <- daily_stats(pm.filled,"PM10",statistic="mean")
## End (Not run)
}
\keyword{data}
\keyword{datagen}
