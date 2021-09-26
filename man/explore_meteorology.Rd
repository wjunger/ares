\name{explore_meteorology}
\alias{explore_meteorology}
\title{Explore Meteorology}
\description{
Explore the exposure-response relationship of meteorology
}
\usage{
explore_meteorology(model, var, df = 4, type = "deviance", new = TRUE, ...)
}
\arguments{
  \item{model}{a model fitted by \code{\link{fit_core}}}
  \item{var}{a quoted string with the name of the meteorology variable}
  \item{df}{degrees of freedom for the smoother. If set to 0 the smooth line is not plotted}
  \item{type}{a quoted string with the type of residuals to plot}
  \item{new}{if \code{TRUE} a new graph window is opened}
  \item{\dots}{{further options for \code{\link{plot}}}}
}
\details{
The scatter plot of residuals against the meteorology variable is smoothed by means of a spline function defined by \code{df}.
}
\value{
This function does not return a value.
}
\references{
Schwartz, J., Spix, C., Touloumi, G. et al. (1996) Methodological issues in studies of air pollution and daily counts of deaths or hospital admissions. \emph{J Epidemiol. Community Health} \bold{50 (suppl 1)}, S12--S18.

Hastie, T., Tibshirani, R. (1990) \emph{Generalized additive models}. 2 ed. Chapman and Hall.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

}

\seealso{\code{\link{plot}},\code{\link{smooth.spline}},\code{\link{setup}}}
\examples{
data(admrio)
setup(admrio,"date")
f <- resp5~s(time)+weekdays+s(humid)
m <- fit_core(f)
explore_meteorology(m,"tmpmax",df=7)
explore_meteorology(m,"humid",df=7)
}
\keyword{ts}
\keyword{smooth}
