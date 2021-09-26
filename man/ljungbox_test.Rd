\name{ljungbox_test}
\alias{ljungbox_test}
\title{Ljung-Box Test}
\description{
Run the Ljung-Box test
}
\usage{
ljungbox_test(x, k = 25, ...)
}
\arguments{
  \item{x}{a vector, usually of residuals}
  \item{k}{an integer indicating the number of lags to compute the autocorrelation coefficients}
  \item{\dots}{further options for \code{\link{acf}}}
}
\details{
This function implements the the Ljung-Box test for autocorrelation coefficients. The null hypothesis of the test is \eqn{H_{0}: \rho(1)=\rho(2)=\cdots=\rho(k)} for \eqn{k<\left(n-1\right)}. The test statistic depends on \code{k}.
}
\value{
A list containing the following objects
  \item{statistic}{the test statistic}
  \item{p.value}{the p-value}
  \item{df}{the degrees of freedom of the chi-squared distribution of the test statistic}
  \item{n.used}{the number of observations used for autocorrelations computation}
  \item{data.name}{the name of the input variable}
}
\references{
Schwartz, J., Spix, C., Touloumi, G. et al. (1996) Methodological issues in studies of air pollution and daily counts of deaths or hospital admissions. \emph{J Epidemiol. Community Health} \bold{50 (suppl 1)}, S12--S18.

Spanos, A. (1999) \emph{Probability Theory and Statistical Inference}. Cambridge.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

}

\seealso{\code{\link{get_residuals}},\code{\link{acf}}}
\examples{
data(admrio)
setup(admrio,"date")
f <- resp5~s(time)+weekdays+s(tmpmax)
m <- fit_core(f)
r <- get_residuals(m)
ljungbox_test(r,25)
}
\keyword{ts}
\keyword{htest}
