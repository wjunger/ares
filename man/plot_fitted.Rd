\name{plot_fitted}
\alias{plot_fitted}
\title{Plot Fitted}
\description{
Plot the observed and predicted values
}
\usage{
plot_fitted(x, gaps = FALSE, date.format = "\%d/\%m/\%Y", new = TRUE, ...)
}
\arguments{
  \item{x}{a model fitted by \code{\link{fit_core}}}
  \item{gaps}{a logical indicating whether \code{.ares.selection} should be honoured}
  \item{date.format}{a string indicating the date format for the horizontal axis. Default is \code{dd/mm/yyyy}. See \code{\link{strptime}} for more options}
  \item{new}{if \code{TRUE} a new graph window is opened}
  \item{\dots}{further options for \code{\link{plot}}}
}
\details{
This function produces a plot with predicted values (points) and fitted ones (line).
}
\value{
This function does not return a value.
}
\references{
McGullagh, P., Nelder, J. A. (1989) \emph{Generalized linear models}. Chapman and Hall.

Hastie, T., Tibshirani, R. (1990) \emph{Generalized additive models}. 2 ed. Chapman and Hall.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

}

\seealso{\code{\link{get_residuals}},\code{\link{gam}},\code{\link{glm}}}
\examples{
data(admrio)
setup(admrio,"date")
f <- resp5~s(time)+weekdays+s(tmpmax)+s(humid)
m <- fit_core(f)
plot_fitted(m)
}
\keyword{ts}
\keyword{regression}
