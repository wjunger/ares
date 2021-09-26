\name{print.summary.pdlm}
\alias{print.summary.pdlm}
\title{Print Summary of PDL Models}
\description{
Method to print a PDLM summary
}
\usage{
\method{print}{summary.pdlm}(x, digits = getOption("digits"), ...)
}
\arguments{
  \item{x}{an object of class \code{summary.pdlm}}
  \item{digits}{an integer indicating the number of decimal places to print. Default is given by the system option \code{digits}}
  \item{\dots}{further option passed by \code{\link{summary.pdlm}}}
}
\value{
This function does not return a value.}
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
summary(dlm)
}
\keyword{ts}
\keyword{regression}
