\name{get_residuals}
\alias{get_residuals}
\title{Get Residuals}
\description{
Extract adjusted residuals from the model
}
\usage{
get_residuals(model, type = "adj_deviance", plot = FALSE, ...)
}
\arguments{
  \item{model}{a model fitted by \code{\link{fit_core}}}
  \item{type}{a quoted string indicating which type of residual to extract. Default is \code{"adj_deviance"}. See Details}
  \item{plot}{a logical indicating whether the residuals should be plotted. See \code{\link{plot_residuals}}}
  \item{\dots}{further options for \code{\link{residuals}}}
}
\details{
The argument \code{type} may be either \code{"deviance"}, \code{"std_deviance"}, \code{"std_scl_deviance"} or \code{"adj_deviance"}. Each of them behaves as described bellow.

\emph{deviance}: Deviance residuals are estimated by \eqn{r_{t}=sign(y_{t}-E(y_{t}))*\sqrt(d_{t})}, where \eqn{d_{t}} is the deviance contribution of the \eqn{t}-th observation. See \code{\link{deviance}} for details on deviance component extraction.

\emph{std\_deviance}: The deviance component is divided by \eqn{(1-h_{t})}, where \eqn{h_{t}} is the \emph{t}-th element of the diagonal of the pseudo hat matrix of the approximating linear model. So they turn into \eqn{r_{t}=sign(y_{t}-E(y_{t}))*\sqrt(d_{t}/(1-h_{t}))}.

\emph{std\_scl\_deviance}: Just like the last one except for the dispersion parameter in its expression, so they have the form \eqn{r_{t}=sign(y_{t}-E(y_{t}))*\sqrt(d_{t}/\phi*(1-h_{t}))}, where \eqn{\phi} is the estimated dispersion parameter of the model. See \code{\link{dispersion}} for \eqn{\phi} estimation.

\emph{adj\_deviance}: These are the deviance residuals multiplied by the estimated coefficient of skewness of the distribution. Thus, for a Poisson model they take the form \eqn{r_{t}=sign(y_{t}-E(y_{t}))*\sqrt(d_{t})*K_{t}}, where \eqn{K_{t}} is given by \eqn{K_{t}=1/(6\sqrt(E(y_{t}))}.

Pierce and Shafer (1986) propose the use of the adjusted deviance residuals over other type of residuals.
}
\value{
A vector of class \code{residuals} with extracted and adjusted residuals of the model.
}
\references{
McGullagh, P., Nelder, J. A. (1989) \emph{Generalized linear models}. Chapman and Hall.

Hastie, T., Tibshirani, R. (1990) \emph{Generalized additive models}. 2 ed. Chapman and Hall.
 
Pierce, D. A., Schafer, D. W. (1986) Residuals in generalized linear models. \emph{Journal of the American Statistical Association}, \bold{81(396)},977--986.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

}

\seealso{\code{\link{residuals}},\code{\link{gam}},\code{\link{glm}}}
\examples{
data(admrio)
setup(admrio,"date")
f <- resp5~s(time)+weekdays+s(tmpmax)+s(humid)
m <- fit_core(f)
r <- get_residuals(m,"adj_deviance")
## plot using appropriate method
plot(r)
}
\keyword{ts}
\keyword{regression}
