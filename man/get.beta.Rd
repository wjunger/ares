\name{get.beta}
\alias{get.beta}
\title{Get Coefficient}
\description{
Transform the constrained coefficients back to the unconstrained form
}
\usage{
get.beta(coeff, var.coeff, lags, degrees, prefix = "beta")
}
\arguments{
  \item{coeff}{a vector of all the coefficients from the distributed lag model}
  \item{var.coeff}{a matrix of variances and covariances of the coefficients from the distributed lag model}
  \item{lags}{an integer indicating the lags used for the distributed lag structure}
  \item{degrees}{an integer indicating the degrees used for the distributed lag structure}
  \item{prefix}{a quoted string to prefix the coefficients names}
}
\details{
This function is used internally by \code{\link{pdlm}} to put the constrained polynomial distributed lag model coefficients, say \eqn{\eta}, back to its original scale \eqn{\beta}.
}
\value{
The function returns a list containing
  \item{beta}{the unconstrained coefficients}
  \item{se.beta}{the unconstrained coefficients standard errors}
}
\references{
Schwartz, J. (2000) The distributed lag between air pollution and daily deaths. \emph{Epidemiology} \bold{11(3)}, 320--326.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

}

\seealso{\code{\link{pdlm}}}

\keyword{ts}
\keyword{regression}
