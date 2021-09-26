\name{print.pdlm}
\alias{print.pdlm}
\title{Print a PDL Model}
\description{
Print some information on a PDL model
}
\usage{
\method{print}{pdlm}(x, digits = getOption("digits"), ...)
}
\arguments{
  \item{x}{a model of class \code{pdlm} created by \code{\link{pdlm}}}
  \item{digits}{an integer indicating the number of decimal places to print. Default is given by the system option \code{digits}}
  \item{\dots}{further options for \code{\link{print}}}
}
\value{
This function does not return a value.
}
\references{
Schwartz, J. (2000) The distributed lag between air pollution and daily deaths. \emph{Epidemiology} \bold{11(3)}, 320--326.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}
}

\seealso{\code{\link{pdlm}}}
\examples{
data(admrio)
setup(admrio,"date")
f <- resp5~s(time)+weekdays+s(tmpmax)+s(humid)
m <- fit_core(f)
dlm <- pdlm(m,"pm10",lags=5,degrees=2)
print(dlm)
}
\keyword{ts}
\keyword{regression}
