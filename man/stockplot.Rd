\name{stockplot}
\alias{stockplot}
\title{Stock Plot}
\description{
Produce a stock type plot
}
\usage{
stockplot(mid, low, high, ref.line = NULL, xlabels = seq(1:length(mid)), ticks = 20,
	mid.pch = 19, lim.pch = 15, graph.scale = NULL, ylim = NULL, ...)
}
\arguments{
  \item{mid}{a vector of mid points for the bars, usually mean values}
  \item{low}{a vector of low points for the bars, usually minimum values}
  \item{high}{a vector of high points for the bars, usually maximum values}
  \item{ref.line}{the position on the y-axis to draw a horizontal line. If is set to \code{NULL} no line is drawn. Default is \code{NULL}}
  \item{xlabels}{a vector of alternate labels for the bars}
  \item{ticks}{an integer indicating the number of tick marks. Default is 20}
  \item{mid.pch}{a number or character indicating the symbol to be used for mid points. Default is a filled circle}
  \item{lim.pch}{a number or character indicating the symbol to be used for end points. Default is a filled square}
  \item{graph.scale}{can be either a logical or a vector with the axis limits. If \code{TRUE} or a vector all the graphs will share the same y-axis scale}
  \item{ylim}{limits for the y-axis scale. If omitted it is computed from values}
  \item{\dots}{further options for \code{\link{plot}}}
}
\details{
This function implements a stock type plot. Often used to represent stock prices closure, it can indicate the average value at the centre of a bar and minimum and maximum values at the ends of the bar. This plot is very useful to represent statistics and their confidence intervals in a convenient way.

This function is used by \code{\link{plot_risk}}.
}
\value{
This function does not return a value.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

}

\seealso{\code{\link{plot}}}
\examples{\dontrun{
data(admrio)
setup(admrio,"date")
m <- c(mean(resp5),mean(resp65))
l <- c(min(resp5),min(resp65))
h <- c(max(resp5),max(resp65))
stockplot(m,l,h,xlabels=c("Resp 5","Resp 65"),xlab="Outcome",ylab="Counts")
}}
\keyword{hplot}
