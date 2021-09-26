\name{unload}
\alias{unload}
\title{Unload}
\description{
Unload the library and do some housekeeping
}
\usage{
unload()
}
\details{
If the library is loaded, this function will detach every instance of \code{.ares.active.dataset} from the search path and remove it from \code{.GlobalEnv}. Then, it will do some memory cleaning and detach every instance of the library.

After \code{unload} is called the library has to be reloaded with \code{\link{library}} .
}
\value{
This function does not return a value.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

}

\seealso{\code{\link{plot}},\code{\link{smooth.spline}},\code{\link{setup}}}
\examples{
## Not run:
# library(ares2)
# unload()
## End (Not run)

}
\keyword{utilities}
