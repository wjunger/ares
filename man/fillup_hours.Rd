\name{fillup_hours}
\alias{fillup_hours}
\title{Fill up Hours}
\description{
Fill up missing hours in the hourly measured datasets
}
\usage{
fillup_hours(input, date.time, readings = 24, cycle = seq(0, 23, by = 24/readings),
	start.date = NULL, end.date = NULL, date.format = "\%d/\%m/\%Y", hour.format =
	"\%H:\%M:\%S", offset = 0, sep = " ")
}
\arguments{
  \item{input}{a data frame. See Details}
  \item{date.time}{a quoted string with the name of a date and time formatted variable}
  \item{readings}{number of measures taken a day. Default is 24 (hourly)}
  \item{cycle}{a vector with the hours measures are taken. Default is a sequence from 0 to 24/\code{measures}}
  \item{start.date}{initial date for the dummy dataset}
  \item{end.date}{final date for the dummy dataset}
  \item{date.format}{date format. See \code{\link{strptime}} for options}
  \item{hour.format}{hour format. See \code{\link{strptime}} for options}
  \item{offset}{an integer with the number of minutes to shift minutes to \code{"00"}. Default is 0}
  \item{sep}{character to separate date and time formats. Default is a space character}
}
\details{
This is an environmental-specific function. Suppose, for instance, that during the day, 24 measurements of a pollutant concentration are taken and if it is missing, the monitor keep the entry out of the log. Thus, at the end of the day less than 24 hourly entries are recorded. However, one needs a regular 24-hour vector to carry on with daily statistics computations, i.e., means, quantiles, etc.

This function reads the input dataset from a data frame. A complete dummy dataset with \code{NA} entries is built from the input parameters. Finally, the original dataset and the dummy dataset are merged using the \code{date.time} variables as reference. The resulting filled up dataset may be output to a data frame or directly to a file which format is supported by \code{\link{export_data}}.
}
\value{
This function returns a data frame with missing filled up with \code{NA}.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

}
\seealso{\code{\link{import_data}},\code{\link{export_data}}}
\examples{
## Not run:
# pm.filled <- fillup_hours(pm_hourly,"Date.Time")
## End (Not run)
}
\keyword{data}
\keyword{datagen}
