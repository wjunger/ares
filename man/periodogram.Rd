\name{periodogram}
\alias{periodogram}
\title{Periodogram}
\description{
Compute and plot the periodogram of the residuals
}
\usage{
periodogram(x, type = "deviance", print = TRUE, rows = 20,
	test = TRUE, new = TRUE, digits = getOption("digits"))
}
\arguments{
  \item{x}{a model fitted by \code{\link{fit_core}} or a vector of residuals}
  \item{type}{a quoted string indicating the type of residuals to extract. Default is \code{"deviance"}. See \code{\link{get_residuals}}}
  \item{print}{a logical indicating whether the statistics should be printed}
  \item{digits}{an integer indicating the number of decimal places to print. Default is given by the system option \code{digits}}
  \item{rows}{an integer giving the number of rows to print}
  \item{test}{a logical indicating whether the uniformity test should be run. See \code{\link{periodogram_test}}}
  \item{new}{if \code{TRUE} a new graph window is opened}
}
\details{
This function computes the Fourier frequencies \eqn{\omega=2\pi j/n} with \eqn{j<n/2} and the periodogram intensity \eqn{I\left(\omega\right)}. A graph of \eqn{I\left(\omega\right)} against the frequencies is plotted.

If \code{object} is a model then the residuals are extracted otherwise the periodogram is computed for the vector supplied.
}
\value{
A data frame with the frequencies, periods, and intensities sorted by the latest in decreasing order.
}
\references{
Schwartz, J., Spix, C., Touloumi, G. et al. (1996) Methodological issues in studies of air pollution and daily counts of deaths or hospital admissions. \emph{J Epidemiol. Community Health} \bold{50 (suppl 1)}, S12--S18.

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
periodogram(m)
}
\keyword{ts}
\keyword{hplot}
\keyword{regression}
