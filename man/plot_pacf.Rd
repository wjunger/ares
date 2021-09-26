\name{plot_pacf}
\alias{plot_pacf}
\title{Plot Partial Autocorrelation Function}
\description{
Plot the autocorrelation function of the residuals
}
\usage{
plot_pacf(x, lags = 25, acf.too = FALSE, type = "deviance", new = TRUE, ...)
}
\arguments{  
  \item{x}{a model fitted by \code{\link{fit_core}}}
  \item{lags}{an integer indicating the number of lags to compute the partial autocorrelation coefficients}
  \item{acf.too}{a logical indicate whether the autocorrelation coefficients should be computed too. Default is \code{FALSE}}
  \item{type}{a quoted string indicating the type of residuals to extract. Default is \code{"deviance"}. See \code{\link{get_residuals}}}
  \item{new}{if \code{TRUE} a new graph window is opened}
  \item{\dots}{further options for \code{\link{pacf}}}
}
\details{
This function computes the partial autocorrelation function (ACF) from the residuals \code{x} and plot it up to \code{lags}. This function was first designed to plot partial autocorrelation coefficients only. However, as some users prefer to analyse de autocorrelation coefficient, the latest was added as option \code{acf.too}.
}
\value{
The function returns a list with the autocorrelations coefficients up to \code{lags}
  \item{acf}{a vector of autocorrelations coefficients. If \code{acf.too} is set to \code{FALSE} it is set to \code{NULL}}
  \item{pacf}{a vector of partial autocorrelations coefficients}
}
\references{
Box, G., Jenkins, G., Reinsel, G. (1994) \emph{Time Series Analysis : Forecasting and Control}. 3 ed., Prentice Hall.

Diggle, P. J. (1990) \emph{Time Series: A Biostatistical Introduction}. Oxford University Press.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

}

\seealso{\code{\link{fit_core}},\code{\link{get_residuals}}}
\examples{
data(admrio)
setup(admrio,"date")
f <- resp5~s(time)+weekdays+s(tmpmax)+s(humid)
m <- fit_core(f)
plot_pacf(m)
}
\keyword{ts}
\keyword{hplot}
\keyword{regression}
