\name{pdlm}
\alias{pdlm}
\title{Polynomial Distributed Lag Models}
\description{
Fit a constrained polynomial distributed lag model
}
\usage{
pdlm(model, var, lags = 5, degrees = 2, ...)
}
\arguments{
  \item{model}{a model fitted by \code{\link{fit_core}}}
  \item{var}{a quoted string with the name with the exposure variable}
  \item{lags}{an integer indicating the number of lags to estimate the effects. Default is 5}
  \item{degrees}{an integer indicating the number of degrees for the constrained polynomial. Default is 2 for a parabolic shape}
  \item{...}{arguments passed on to other methods. See \code{\link{family}}}
}
\details{
This function updates \code{model} with the unconstrained distributed lag models using \code{\link{pdl}}. Then, the unconstrained coefficients and their standard errors are extracted using \code{\link{get.beta}}.

This model is thoroughly discussed in Schwartz (2000).
}
\value{
The class \code{pdlm} is added to the model inheritance and the following list is returned
	\item{cmodel}{the fitted constrained model}
	\item{variate}{the vector with exposure variate data}
	\item{var.name}{the name of exposure variate}
	\item{beta}{the unconstrained coefficients}
	\item{lags}{an integer indicating the lags used for the distributed lag structure}
	\item{degrees}{an integer indicating the degrees used for the polynomial in the distributed lag structure}
	\item{call}{function call}
}
\references{
Schwartz, J. (2000) The distributed lag between air pollution and daily deaths. \emph{Epidemiology} \bold{11(3)}, 320--326.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

}

\seealso{\code{\link{gam}},\code{\link{glm}},\code{\link{fit_core}}}
\examples{
data(admrio)
setup(admrio,"date")
f <- resp5~s(time)+weekdays+s(tmpmax)+s(humid)
m <- fit_core(f)
dlm <- pdlm(m,"pm10",lags=5,degrees=2)
}
\keyword{ts}
\keyword{regression}
