\name{dispersion}
\alias{dispersion}
\title{Dispersion Parameter}
\description{
Estimate dispersion parameter
}
\usage{
dispersion(model)
}
\arguments{
  \item{model}{a model fitted by \code{\link{fit_core}}}
}
\details{
This is an consistent and asymptotic estimation of the dispersion parameter based on generalized Pearson statistics for GLM models. It is given \eqn{\frac{\sum_{t=1}^{n}\left(y_{t}-E\left(y_{t}\right)\right)^{2}/E\left(y_{t}\right)}{n-df}}, where \eqn{n} is the number of observations and \eqn{df} is the number of degrees of freedom used by the model.
}
\value{
A scalar with the statistic.
}
\references{
McGullagh, P., Nelder, J. A. (1989) \emph{Generalized linear models}. Chapman and Hall.

Hastie, T., Tibshirani, R. (1990) \emph{Generalized additive models}. 2 ed. Chapman and Hall.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

}

\seealso{\code{\link{fit_core}},\code{\link{gam}},\code{\link{glm}}}
\examples{
data(admrio)
setup(admrio,"date")
f <- resp5~s(time)+weekdays+s(tmpmax)+s(humid)
m <- fit_core(f)
dispersion(m)
}
\keyword{regression}
