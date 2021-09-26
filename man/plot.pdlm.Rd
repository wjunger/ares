\name{plot.pdlm}
\alias{plot.pdlm}
\title{Plot of PDL Models}
\description{
Plot the effects estimated with a polynomial distributed lag model
}
\usage{
\method{plot}{pdlm}(x, unit = 10, confidence.level = 0.95,
	labels = NULL, new = TRUE, ...)
}
\arguments{
  \item{x}{a model of class \code{pdlm} created by \code{\link{pdlm}}}
  \item{unit}{a vector indicating the units for relative risk computation. Default is 10. See \code{\link{estimate_risks}} for more options}
  \item{confidence.level}{confidence level for interval computation}
  \item{labels}{a vector of quoted strings with alternate labels for the pollutants. Default is the time lag index. See Details}
  \item{new}{if \code{TRUE} a new graph window is opened}
  \item{\dots}{further options for \code{\link{stockplot}}}
}
\details{
This function extracts the coefficients and their standard errors from a polynomial distributed lag model in \code{x} and plots a proper graph.

If \code{labels} is a vector, it must have the same length as the number of lags plus 1; if it is a quoted string it will be used as a prefix for the labels with the time lag index; and if it is \code{NULL} automatically created labels \dQuote{Lag i} will be used.
}
\value{
This function does not return a value.
}
\references{
Schwartz, J. (2000) The distributed lag between air pollution and daily deaths. \emph{Epidemiology} \bold{11(3)}, 320--326.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}
}

\seealso{\code{\link{pdlm}},\code{\link{estimate_risks}},\code{\link{plot_risk}}}
\examples{
data(admrio)
setup(admrio,"date")
f <- resp5~s(time)+weekdays+s(tmpmax)+s(humid)
m <- fit_core(f)
dlm <- pdlm(m,"pm10",lags=5,degrees=2)
plot(dlm)
}
\keyword{ts}
\keyword{regression}
