\name{sincos}
\alias{sincos}
\title{Sinusoidal Basis Generating Function}
\description{
Generate a basis for a sinusoidal function at some fixed period
}
\usage{
sincos(period, n = NULL, largest.period = 365)
}
\arguments{
  \item{period}{a numeric value indicating the seasonal pattern period}
  \item{n}{an integer indicating the number of observations}
  \item{largest.period}{a numeric value indicating the largest seasonal pattern period. Default is an year or 365 days}
}
\details{
This function returns a two-column matrix containing both sine and cosine component evaluated at the given period. The columns are named accordingly.

The function is designed to be used directly from within a formula. Therefore, there is no need to create the sine and cosine vectors prior to the fitting process.

When using with the \bold{ares} library, \code{n} may be omitted and it will be taken from the length of \code{ares.selection}. 
}
\value{
A two-column matrix containing sine and cosine
}
\references{
Schwartz, J., Spix, C., Touloumi, G. et al. (1996) Methodological issues in studies of air pollution and daily counts of deaths or hospital admissions. \emph{J Epidemiol. Community Health} \bold{50 (suppl 1)}, S12--S18.

Box, G., Jenkins, G., Reinsel, G. (1994) \emph{Time Series Analysis : Forecasting and Control}. 3 ed., Prentice Hall.

Diggle, P. J. (1990) \emph{Time Series : A Biostatistical Introduction}. Oxford University Press.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

}

\seealso{\code{\link{fit_core}},\code{\link{get_residuals}}}
\examples{
data(admrio)
setup(admrio,"date")
f <- resp5~s(time)+weekdays+s(tmpmax)+s(humid)+sincos(90)+sincos(180)
m <- fit_core(f)
print_summary(m)
}
\keyword{ts}
\keyword{hplot}
\keyword{regression}
