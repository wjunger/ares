\name{lagvar}
\alias{lagvar}
\title{Matrix of lagged variable}
\description{
Create a matrix of a variable lagged up to some instant
}
\usage{
lagvar(var, k)
}
\arguments{
  \item{var}{a vector of the variable to lag}
  \item{k}{an integer indicating the number of lags to compute}
}
\details{
This function creates a matrix with \eqn{k+1} columns. The columns hold \code{var} lagged from 0 to \eqn{k}.

It is designed to be called from within \code{\link{pdl}}.
}
\value{
A matrix with the lagged variate.
}
\references{
Schwartz, J. (2000) The distributed lag between air pollution and daily deaths. \emph{Epidemiology} \bold{11(3)}, 320--326.

Diggle, P. J. (1990) \emph{Time Series: A Biostatistical Introduction}. Oxford University Press.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

}

\seealso{\code{\link{pdlm}}}

\keyword{ts}
\keyword{regression}
