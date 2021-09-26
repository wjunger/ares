\name{diagnostics}
\alias{diagnostics}
\title{Diagnostics}
\description{
Perform a series of model dignostics
}
\usage{
diagnostics(model, single.graph = TRUE)
}
\arguments{
  \item{model}{a model fitted by \code{\link{fit_core}}}
  \item{single.graph}{a logical indicating whether all graphs should be plotted on the same device windows}
}
\details{
This function is a short-cut to several common diagnostics used when analysing time series of outcomes in environmental epidemiology. The following procedures are executed by \code{diagnostics}: plot of fitted values with \code{\link{plot_fitted}}, plot of residuals with \code{\link{plot_residuals}}, plot of Cook's distance with \code{\link{plot_cook}}, plot of partial autocorrelation function with \code{\link{plot_pacf}}, plot of periodogram with \code{\link{periodogram}}, quantile-quantile plot with \code{\link{plot_qq}}, and information on parameters and smooth functions estimates \code{\link{print_summary}}.

If all the graphs are plotted on the same device window, it can be rather difficult to read, but it is more suitable for publications and the graph can be saved with high quality format using \code{\link{save.plot}}.
}
\value{
This function does not return a value.
}
\references{
McGullagh, P., Nelder, J. A. (1989) \emph{Generalized linear models}. Chapman and Hall.

Hastie, T., Tibshirani, R. (1990) \emph{Generalized additive models}. 2 ed. Chapman and Hall.
 
Pierce, D. A., Schafer, D. W. (1986) Residuals in generalized linear models. \emph{Journal of the American Statistical Association}, \bold{81(396)},977--986.

Schwartz, J., Spix, C., Touloumi, G. et al. (1996) Methodological issues in studies of air pollution and daily counts of deaths or hospital admissions. \emph{J Epidemiol. Community Health} \bold{50 (suppl 1)}, S12--S18.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

}

\seealso{\code{\link{fit_core}},\code{\link{summary}}}
\examples{
data(admrio)
setup(admrio,"date")
f <- resp5~s(time)+weekdays+s(tmpmax)+s(humid)
m <- fit_core(f)
diagnostics(m)
}
\keyword{hplot}
\keyword{regression}
