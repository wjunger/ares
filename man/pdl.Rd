\name{pdl}
\alias{pdl}
\title{Polynomial Distributed Lag Basis}
\description{
Set up a basis for a constrained polynomial distributed lag model
}
\usage{
pdl(var, lags, degrees)
}
\arguments{
  \item{var}{a vector with the variate to constrain on the polynomial functional}
  \item{lags}{an integer indicating the maximum lags for the basis}
  \item{degrees}{an integer indicating the degrees used for the polynomial in the distributed lag structure}
}
\details{
Each column of the resulting matrix but the first one will be the form \eqn{W_d=Z_1+2^dZ^2+..+q^dZ_q} and  the first one \eqn{W_0=Z_0 + Z_1+..+ Z_q}, where \eqn{Z_i} is \code{var} at lag \eqn{i}.

This function is intended to be called internally by \code{\link{pdlm}}.
}
\value{
This function returns a matrix with \code{degrees+1} columns.
}
\references{
Schwartz, J. (2000) The distributed lag between air pollution and daily deaths. \emph{Epidemiology} \bold{11(3)}, 320--326.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

}

\seealso{\code{\link{pdlm}}}

\keyword{ts}
\keyword{regression}
