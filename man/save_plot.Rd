\name{save_plot}
\alias{save_plot}
\title{Save Plot}
\description{
Conveniently save the active plot window to some known graph file format
}
\usage{
save_plot(file, width = 520, height = 480)
}
\arguments{
  \item{file}{quoted string with the file name}
  \item{width}{width in pixels}
  \item{height}{height in pixels}
}
\details{
This function tries to guess the filter to use from the file name extension. Supported graph formats are portable bitmap (png), jpeg (jpg), portable document format (pdf), encapsulated post script (eps), and scalable vector graphics (svg). If \code{off} is set to \code{TRUE}, the active window device is closed after the graph is saved.

Alternatively, one can use the \R native \code{\link{savePlot}} that essentially does the same job.

}
\value{ 
This function does not return a value.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

}
\seealso{\code{\link{dev.copy}},\code{\link{dev.copy2eps}}}
\examples{
data(admrio)
setup(admrio,"date")
plot_event("resp5")
save_plot("plot_resp5.png")
}
\keyword{device}
