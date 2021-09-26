\name{print_risk}
\alias{print_risk}
\title{Print Risk}
\description{
Print a risk object nicely formatted
}
\usage{
print_risk(x, digits = getOption("digits"), ...)
}
\arguments{
  \item{x}{an object of class \code{risk} or \code{pdlm.risk} output by \code{\link{estimate_risks}}}
  \item{digits}{an integer indicating the number of decimal places to print. Default is given by the system option \code{digits}}
  \item{\dots}{further options to pass to \code{\link{print}}}
}
\details{
This function prints out the risk tables nicely formatted.
}
\value{
This function does not return a value.
}
\references{
Schwartz, J., Spix, C., Touloumi, G. et al. (1996) Methodological issues in studies of air pollution and daily counts of deaths or hospital admissions. \emph{J Epidemiol. Community Health} \bold{50 (suppl 1)}, S12--S18.

McGullagh, P., Nelder, J. A. (1989) \emph{Generalized linear models}. Chapman and Hall.

Hastie, T., Tibshirani, R. (1990) \emph{Generalized additive models}. 2 ed. Chapman and Hall.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

}

\seealso{\code{\link{print}},\code{\link{estimate_risks}}}
\examples{
data(admrio)
setup(admrio,"date")
f <- resp5~s(time)+weekdays+s(tmpmax)+s(humid)
m <- fit_core(f)
rr <- estimate_risks(m,c("pm10","so2"),digits=3,labels=c("PM10","SO2"),
	 method="singlelag",lag.struc=list(l=0:2,ma=1:5))

print(rr)
}
\keyword{regression}
