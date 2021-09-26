\name{lspline}
\alias{lspline}
\title{Linear Spline}
\description{
Generate a basis for a linear spline for piecewise regression
}
\usage{
lspline(var, knots = NULL, nknots = NULL, percentiles = NULL, 
	marginal = FALSE, names = NULL)
}
\arguments{
  \item{var}{the covariate to generate the basis}
  \item{knots}{a vector indicating the positions of the knots. See Details}
  \item{nknots}{an integer indicating the number of knots. See Details}
  \item{percentiles}{a vector of percentiles indicating the positions of the knots. See Details}
  \item{marginal}{logical indicating whether the effect for each interval of the spline represents the change in slope from the preceding one. Default is \code{FALSE}, so the coefficients represent the slope of the interval}
  \item{names}{a vector of quoted strings with alternate names for the columns of the resulting matrix. The length must be the number of knots plus 1}
}
\details{
This functions implements in \R the Stata \code{mkspline} command. For further references see the article Gould (1993).

It is necessary to set either of \code{knots}, \code{nknots}, or \code{percentiles}. If more than one is specified, it will throw an error. If \code{nknots} is specified the intervals created will be equally spaced.

The function is designed to be called directly from within a formula. Therefore, there is no need to create the basis prior to fitting the model.

This function is very useful for exploring the exposure-response curve of the pollutants.
}
\value{
This function returns a matrix with the basis for a linear spline. The number of columns in the matrix is the number of knots plus 1. The columns are named accordingly or after \code{names}
}
\references{
Gould, W. W. (1993) Linear splines and piecewise linear functions. \emph{Stata Technical Bulletin} \bold{15}, 13--17.

Schwartz, J., Spix, C., Touloumi, G. et al. (1996) Methodological issues in studies of air pollution and daily counts of deaths or hospital admissions. \emph{J Epidemiol. Community Health} \bold{50 (suppl 1)}, S12--S18.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

}

\seealso{\code{\link{smooth.spline}},\code{\link{setup}}}
\examples{
data(admrio)
setup(admrio,"date")
%sp.pm10.1 <- lspline(pm10,knots=c(52.50,61.02,74.75))
%sp.pm10.2 <- lspline(pm10,nknots=3)
%sp.pm10.3 <- lspline(pm10,percentiles=c(0.25,0.50,0.75))
%
%## within a formula
f <- resp5~s(time)+weekdays+s(tmpmax)+s(humid)+lspline(pm10,knots=c(52.50,61.02,74.75))
m <- fit_core(f)
print_summary(m)
}
\keyword{regression}
