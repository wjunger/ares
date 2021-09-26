\name{ma}
\alias{ma}
\title{Moving Average}
\description{
Compute simple moving average of a variable.
}
\usage{
ma(var, begin, end, selection = TRUE)
}
\arguments{
  \item{var}{variable to average}
  \item{begin}{initial index}
  \item{end}{final index}
  \item{selection}{a logical indicating whether \code{.ares.selection} should honoured}
}
\details{
This function is used mainly to compute accumulated exposure to meteorological factors and/or pollutants. The arguments \code{begin} and \code{end} control the window size for simple moving average computation, for instance, if one sets this parameters to 0 and 2 respectively, the resulting averages will be based on the concurrent day and the two previous days.

The function may be used inside a model formula to estimate the accumulated exposure effects.
}
\value{
Return \code{var} averaged and if \code{.ares.selection} is set, averaged \code{var} is cut accordingly.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

}

\seealso{\code{\link{mean}},\code{\link{setup}}}
%\examples{
%data(admrio)
%setup(admrio,"date")
%# compute moving average of concurrent and previous two days
%pm10.ma02 <- ma(pm10,0,2)
%}
\keyword{ts}
\keyword{univar}
