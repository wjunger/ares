\name{summary.pdlm}
\alias{summary.pdlm}
\title{Summary of PDL Models}
\description{
Output the summary of a distributed lag model
}
\usage{
\method{summary}{pdlm}(object, ...)
}
\arguments{
  \item{object}{a model of class \code{pdlm} created by \code{\link{pdlm}}}
  \item{\dots}{further options for \code{\link{print.summary.pdlm}}}
}
\details{
The function prints information about the model.
}
\value{
Invisibly returns a list containing
  \item{call}{function call}
  \item{coef.table}{a table containing coefficients, standard errors, t values and p-values}
  \item{resid.stats}{some statistics on the deviance residuals of the constrained model}
  \item{deviance}{the residual deviance  of the constrained model}
  \item{null.deviance}{the null deviance of the constrained model}
  \item{aic}{Akaike Information Criterion of the constrained model}
  \item{iter}{number of Fisher scoring iterations of the constrained model}
  \item{df.null}{null degrees of freedom of the constrained model}
  \item{df residual}{residual degrees of freedom of the constrained model}
  \item{na}{number of observations lost due missing data}
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
summary(dlm)
}
\keyword{ts}
\keyword{regression}
