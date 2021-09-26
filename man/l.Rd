\name{l}
\alias{l}
\title{Lagged Variables}
\description{
Compute lagged variables
}
\usage{
l(var, k, selection = TRUE)
}
\arguments{
  \item{var}{variable for lagging}
  \item{k}{lag}
  \item{selection}{a logical indicating whether \code{.ares.selection} should honoured}
}
\details{
This function is used mainly to assess lagged exposure to meteorological factors and/or pollutants. The argument \code{k} controls the lag distance.

The function may be used inside a model formula to estimate the lagged exposure effects.
}
\value{
Return the variable \code{var} lagged. If \code{.ares.selection} is set, lagged \code{var} is cut accordingly.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

}

\seealso{\code{\link{mean}},\code{\link{setup}}}
%\examples{
%data(admrio)
%setup(admrio,"date")
%# compute lag2 of pm10
%pm10.l2 <- l(pm10,2)
%}
\keyword{ts}
\keyword{univar}
