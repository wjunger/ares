\name{rr.eval}
\alias{rr.eval}
\title{Relative Risk Evaluation}
\description{
Compute relative risk for some unit and confidence intervals
}
\usage{
rr.eval(beta, se, unit = 1, confidence.level = 0.95)
}
\arguments{
  \item{beta}{regression coefficient of the exposure variate}
  \item{se}{regression standard error of the exposure variate}
  \item{unit}{a value by which the coefficient is multiplied}
  \item{confidence.level}{confidence level for interval computation}
}
\details{
The function simply compute the relative risk and the confidence interval for a given amount of pollutant in \code{unit}. It is intended to be called internally by \code{\link{estimate_risks}}.
}
\value{
A vector containing the relative risk and confidence limits.
}
\references{
Schwartz, J., Spix, C., Touloumi, G. et al. (1996) Methodological issues in studies of air pollution and daily counts of deaths or hospital admissions. \emph{J Epidemiol. Community Health} \bold{50 (suppl 1)}, S12--S18.

McGullagh, P., Nelder, J. A. (1989) \emph{Generalized linear models}. Chapman and Hall.

Hastie, T., Tibshirani, R. (1990) \emph{Generalized additive models}. 2 ed. Chapman and Hall.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

}

\seealso{\code{\link{fit_core}},\code{\link{plot_risk}}}
\examples{
data(admrio)
setup(admrio,"date")
f <- resp5~s(time)+weekdays+s(tmpmax)+s(humid)
m <- fit_core(f)
rr <- estimate_risks(m,c("pm10","so2"),digits=3,labels=c("PM10","SO2"),
	method="singlelag",lag.struc=list(l=0:2,ma=1:5))
}
\keyword{models}
\keyword{regression}
