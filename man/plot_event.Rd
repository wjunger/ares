\name{plot_event}
\alias{plot_event}
\title{Plot Event}
\description{
Plot the daily series of event counts.
}
\usage{
plot_event(x, df = 4, gaps = FALSE, type = "p",
	title = NULL, date.format = "\%d/\%m/\%Y", new = TRUE, ...)
}
\arguments{
  \item{x}{a quoted string with the name of the counts variable}
  \item{df}{an integer with the degrees of freedom of the spline. If set to 0, no line is plotted}
  \item{gaps}{a logical indicating whether \code{.ares.selection} should be honoured}
  \item{type}{a string indicating the type of graphic. Use \code{"p"} for points and \code{"l"} for line. Default is \code{"p"}. See \code{\link{plot}}}
  \item{title}{a string supplying a title for the graphic. If \code{NULL}, an automatic title is generated}
  \item{date.format}{a string indicating the date format for the horizontal axis. Default is \code{dd/mm/yyyy}. See \code{\link{strptime}} for more options}
  \item{new}{if \code{TRUE} a new graph window is opened}
  \item{\dots}{further options for \code{\link{plot}}}
}
\details{
This function and \code{\link{plot_pollutant}} are essentially the same, only some default options differ. There are two of them for convenience only.
}
\value{
This function does not return a value.
}
\references{
Schwartz, J., Spix, C., Touloumi, G. et al. (1996) Methodological issues in studies of air pollution and daily counts of deaths or hospital admissions. \emph{J Epidemiol. Community Health} \bold{50 (suppl 1)}, S12--S18.

Hastie, T., Tibshirani, R. (1990) \emph{Generalized additive models}. 2 ed. Chapman and Hall.

Peng, R., Dominici, F. (2008) Statistical methods for environmental epidemiology with R.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

}

\seealso{\code{\link{plot}},\code{\link{smooth.spline}},\code{\link{setup}}}
\examples{
data(admrio)
setup(admrio,"date")
plot_event("resp5",df=7)
}
\keyword{hplot}
\keyword{smooth}
