\name{plot_risk}
\alias{plot_risk}
\title{Plot Risk}
\description{
Plot estimated risks
}
\usage{
plot_risk(x, labels = rownames(x), new = TRUE, graph.scale = FALSE, ...)
}
\arguments{
  \item{x}{an object of class \code{risk} or \code{pdlm.risk} output by \code{\link{estimate_risks}}}
  \item{labels}{a vector of alternate labels for exposures with the same dimension of \code{risks}. Default is the vector of rows names of \code{x}}
  \item{new}{if \code{TRUE}, a new graph window is opened}
  \item{graph.scale}{can be either a logical or a vector with the axis limits. If \code{TRUE} or a vector all the graphs will share the same y-axis scale}
  \item{\dots}{further options for \code{\link{stockplot}}}
}
\details{
This function will plot relative risks computed by \code{\link{estimate_risks}}. It is seldom called directly.
}
\value{
This function does not return a value.
}
\references{
Schwartz, J., Spix, C., Touloumi, G. et al. (1996) Methodological issues in studies of air pollution and daily counts of deaths or hospital admissions. \emph{J Epidemiol. Community Health} \bold{50 (suppl 1)}, S12--S18.

McGullagh, P., Nelder, J. A. (1989) \emph{Generalized linear models}. Chapman and Hall.

Hastie, T., Tibshirani, R. (1990) \emph{Generalized additive models}. 2 ed. Chapman and Hall.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

}

\seealso{\code{\link{plot}},\code{\link{estimate_risks}}}
\examples{
data(admrio)
setup(admrio,"date")
f <- resp5~s(time)+weekdays+s(tmpmax)+s(humid)
m <- fit_core(f)
rr <- estimate_risks(m,c("pm10","so2"),digits=3,labels=c("PM10","SO2"),
	method="singlelag",lag.struc=list(l=0:2,ma=1:5))

plot(rr)
}
\keyword{hplot}
\keyword{regression}
