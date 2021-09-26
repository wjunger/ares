\name{plot_residuals}
\alias{plot_residuals}
\title{Plot Residuals}
\description{
Extract and/or plot the model residuals
}
\usage{
plot_residuals(x, gaps = FALSE, type = "deviance",
	 band=c(-3,3), date.format = "\%d/\%m/\%Y", new = TRUE, ...)
}
\arguments{
  \item{x}{a model fitted by \code{\link{fit_core}} or a vector of the class \code{residuals}}
  \item{gaps}{a logical indicating whether \code{.ares.selection} should be honoured}
  \item{type}{a quoted string indicating the type of residuals to extract. Default is \code{"deviance"}. See \code{\link{get_residuals}}}
  \item{band}{a vector limits to plot the band. If \code{NULL}, the band is omitted. Default is \code{c(-3,3)}}
  \item{date.format}{a string indicating the date format for the horizontal axis. Default is \code{dd/mm/yyyy} See \code{\link{strptime}} for more options}
  \item{new}{if \code{TRUE} a new graph window is opened}
  \item{\dots}{further options for \code{\link{plot}}}
}
\details{
If \code{x} is a model then the residuals are extracted with \code{\link{get_residuals}}. If \code{x} is a vector of the class \code{residuals} the function works as a method for the generic function \code{\link{plot}}, therefore only the graph is plotted.
}
\value{
A vector of class \code{residuals} with extracted and/or adjusted residuals of the model.
}
\references{
McGullagh, P., Nelder, J. A. (1989) \emph{Generalized linear models}. Chapman and Hall.

Hastie, T., Tibshirani, R. (1990) \emph{Generalized additive models}. 2 ed. Chapman and Hall.
 
Pierce, D. A., Schafer, D. W. (1986) Residuals in generalized linear models. \emph{Journal of the American Statistical Association}, \bold{81(396)},977--986.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}
}

\seealso{\code{\link{get_residuals}},\code{\link{resid}},\code{\link{gam}},\code{\link{glm}}}
\examples{
data(admrio)
setup(admrio,"date")
f <- resp5~s(time)+weekdays+s(tmpmax)+s(humid)
m <- fit_core(f)
## extract the residuals and plot
r1 <- plot_residuals(m)

## extract first then plot
r <- get_residuals(m,"adj_deviance")
plot(r)
}
\keyword{hplot}
\keyword{regression}
