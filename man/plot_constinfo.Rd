\name{plot_constinfo}
\alias{plot_constinfo}
\title{Plot Constant Information}
\description{
Plot constant information on the response scale
}
\usage{
plot_constinfo(x, type = "deviance", new = TRUE, ...)
}
\arguments{
  \item{x}{a model fitted by \code{\link{fit_core}}}
  \item{type}{a quoted string indicating the type of residuals to extract. Default is \code{"deviance"}. See \code{\link{get_residuals}}}
  \item{new}{if \code{TRUE} a new graph window is opened}
  \item{\dots}{further options for \code{\link{plot}}}
}
\details{
This function extracts the residuals of \code{x} and generates a plot of the residuals against 2 times the square root of the fitted values.
}
\value{
This function does not resturn a value.
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
plot_constinfo(m)
}
\keyword{hplot}
\keyword{regression}
