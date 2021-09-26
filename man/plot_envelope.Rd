\name{plot_envelope}
\alias{plot_envelope}
\title{Plot Envelope}
\description{
Plot a simulated envelope of the residuals}
\usage{
plot_envelope(x, rep = 39, type = "deviance", new = TRUE, ...)
}
\arguments{
  \item{x}{a model fitted by \code{\link{fit_core}}}
  \item{rep}{an integer indicating the number of replications. Default is 19. See Details}
  \item{type}{a quoted string indicating the type of residuals to extract. Default is \code{"deviance"}. See \code{\link{get_residuals}}}
  \item{new}{if \code{TRUE} a new graph window is opened}
  \item{\dots}{further options for \code{\link{qqnorm}}}
}
\details{
This function implements a simulated envelope of residuals. It is a plot to check if the normality assumptions of the residuals hold. Although the authors originally suggested 39 replications, at least a hundred may be reasonable, since the processing power of modern computers is no longer a limitation. However, the re-sampling strategy of the method is still time consuming.

The algorithm is implemented for Poisson family models only. Although it supports \emph{quasi-poisson} models, the results may be inaccurate.
}
\value{
This function does not return a value.
}
\references{
McGullagh, P., Nelder, J. A. (1989) \emph{Generalized linear models}. Chapman and Hall.

Hastie, T., Tibshirani, R. (1990) \emph{Generalized additive models}. 2 ed. Chapman and Hall.
 
Cook, R. D., Weisberg, S. (1982) \emph{Residuals and Influence in Regression}. Chapman and Hall.

Atkinson, A.C. (1985) \emph{Plots, Transformations and Regression}. Oxford University Press.

Pierce, D. A., Schafer, D. W. (1986) Residuals in generalized linear models. \emph{Journal of the American Statistical Association}, \bold{81(396)},977--986.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}
}

\seealso{\code{\link{fit_core}},\code{\link{get_residuals}}}
\examples{
data(admrio)
setup(admrio,"date")
f <- resp5~s(time)+weekdays+s(tmpmax)+s(humid)
m <- fit_core(f)
plot_envelope(m)
}
\keyword{hplot}
\keyword{regression}
