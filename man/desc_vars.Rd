\name{desc_vars}
\alias{desc_vars}
\title{Describe Variables}
\description{
Compute some common descriptives statistics of variables
}
\usage{
desc_vars(vars, by = NULL, stats = c("n", "na", "mean", "sd", 
	"min", "max", "centiles"), probs = c(0.25, 0.5, 0.75), 
	labels = NULL, print = TRUE, digits = getOption("digits"), ...)
}
\arguments{
  \item{vars}{a vector of quoted strings with the variables names or a formula, e.g. \code{~x+y}}
  \item{by}{a factor variable for subset analysis.It will be coerced to a factor. Default is \code{NULL}}
  \item{stats}{a vector of quoted strings with the statistics to be computed. See Details}
  \item{probs}{a vector of probabilities for the quantiles computation}
  \item{labels}{a vector of quoted strings with alternate labels for the variables. If omitted \code{vars} is used}
  \item{print}{a logical indicating whether the statistics should be printed}
  \item{digits}{an integer indicating the number of decimal places to print. Default is given by the system option \code{digits}}
  \item{\dots}{further arguments for \code{\link{mean}}, basically, support for trimming}
}
\details{
By default, this function computes the following statistics: \emph{n}: number of observations, \emph{na}: number of missing values, \emph{mean}: variable means, \emph{sd}: standard deviation, \emph{min}: minimum value; \emph{max}: maximum value, and \emph{centiles}: quantiles controlled by the vector \code{probs}. Any subset of these statistics can be selected. However, if other than these is indicated it will throw an error. The statistics in \code{stats} are not case sensitive.

If the argument \code{by} is set, the statistics are computed for each subset defined by the supplied vector.

These statistics are frequently found on publications regarding air pollution and health effects.
}
\value{                                      
The function invisibly returns a data frame with a variable per line and a statistic per column.
}
\references{
Schwartz, J., Spix, C., Touloumi, G. et al. (1996) Methodological issues in studies of air pollution and daily counts of deaths or hospital admissions. \emph{J Epidemiol. Community Health} \bold{50 (suppl 1)}, S12--S18.

Whetherill, G. B. \emph{Elementary Statistical Methods}. 2 ed. Chapman and Hall.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

}
\seealso{\code{\link{summary}},\code{\link{mean}},\code{\link{quantile}}}
\examples{
data(admrio)
setup(admrio,"date")
## no stratified analysis
d <- desc_vars(~resp5+pm10+so2+co,stats=c("n","mean","sd","centiles"),digits=3)
## stratified analysis
d <- desc_vars(~resp5+pm10+so2+co,by=quarters,stats=c("n","mean","sd","centiles"),
	digits=3)
}
\keyword{ts}
\keyword{univar}
