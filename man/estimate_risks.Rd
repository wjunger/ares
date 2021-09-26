\name{estimate_risks}
\alias{estimate_risks}
\title{Estimate Risks}
\description{
Estimate the effects of the pollutants
}
\usage{
estimate_risks(model,pollutant,unit=10,confidence.level=.95,method="singlelag",
	perc.rr=TRUE,interaction=NULL,lag.struc=list(l=0:5,ma=NULL),
	pdlm.struc=list(l=5,d=2,overall=TRUE),overdispersion=FALSE,labels=NULL,
	print=TRUE,digits=getOption("digits"),plot=TRUE,new=TRUE,graph.scale=FALSE,
	verbose=TRUE,...)
}
\arguments{
  \item{model}{a model fitted by \code{\link{fit_core}}}
  \item{pollutant}{a vector with the names of the variables to estimate the effects or a formula, e.g. \code{~x+y}}
  \item{unit}{a single value or a vector indicating the units for relative risk computation. Default is 10 for all pollutants. See Details}
  \item{confidence.level}{confidence level for interval computation}
  \item{method}{a string indicating the method for effect estimation. Default is \code{"singlelag"}. See Details}
  \item{perc.rr}{logical. If \code{TRUE} the effects are reported in terms of percent changes in risk. Default is \code{TRUE}}
  \item{interaction}{a string indicating the 2-level interaction term usually used for seasonal effects estimation. Default is \code{NULL}. See Details}
  \item{lag.struc}{a list with the single lag model structure. Default is lag up to 5 days. See Details}
  \item{pdlm.struc}{a list with the polynomial distributed lag model structure. Default is lag up to 5 and a 2-degree polynomial function. If \code{overall=TRUE} overall effect is computed. See Details}
  \item{overdispersion}{a logical indicating whether confidence intervals should account for the extra variability}
  \item{labels}{a vector of quoted strings with alternate labels for the pollutants. Default is the names of the variables in \code{pollutant}}
  \item{print}{a logical indicating whether the statistics should be printed}
  \item{digits}{an integer indicating the number of decimal places to print. Default is given by the system option \code{digits}}
  \item{plot}{a logical indicating whether the estimated risks should be plotted. Default is \code{TRUE}. See \code{\link{plot_risk}}}
  \item{new}{if \code{TRUE} a new graph window is opened}
  \item{graph.scale}{can be either a logical or a vector with the axis limits. If \code{TRUE} or a vector all the graphs will share the same y-axis scale}
  \item{verbose}{a logical indicating whether extra information should be printed during the iterations}
  \item{\dots}{further options for \code{\link{plot_risk}}}
}
\details{
This function estimates the effects for each pollutant in \code{pollutant} using the estimation approach set in \code{method}. If \code{method} is set to either \code{"singlelag"} or \code{"dual"} the effects are estimated independently for each exposure. If it is set to \code{"pdlm"} then a polynomial distributed lag model is used to estimate the effects using the lag structure passed to \code{\link{pdlm}} using \code{pdlm.struc} option.

The \code{lag.struc} argument is a list containing the number of lags (l) and/or the moving averages (ma) of the pollutants as the measure of exposure. This mode allows more than one pollutant at a time. The general list structure is \code{list(l=,ma=,ma.base=,labels=)}, where \code{l} and \code{ma} are vectors indicating the lagged exposures in order to estimate the effects. If \code{ma.base} is omitted or \code{ma.base=0}, the moving averages will range from the concurrent day to each element in \code{ma}. If \code{labels} is missing or set to \code{NULL}, a generic label will be used.

The \code{pdlm.struc} argument is a list containing the number of lags (l) and the degrees (deg) to be passed to \code{\link{pdlm}}. This mode allows more than one pollutant at a time. The general list structure is \code{list(l=,deg=,labels=)}, where \code{l} is an integer and \code{deg} is an integer or vector with the same length as \code{pollutant}. If \code{deg} is an integer all the pollutants will share the same polynomial structure. If \code{labels} is missing or set to \code{NULL}, a generic label will be used instead. Overall estimate and confidence interval will be plotted if \code{overall=TRUE}. 

Dual pollutant models may be estimated by setting the option \code{method} to \code{"dual"}. This method will estimate a model for each combination of two of the pollutants set in \code{pollutant}. One can set the option \code{lag.struc} the same way as in the single lag models, however both pollutants in the model will share the same lag structure, i.e., the effects of both pollutant will be estimated by using the same lagged exposure. It is not a serious limitation though. Due to the manner the pollutants effects matrices are stored in the risk array, one should avoid reading them directly. Use the \code{\link{print_risk}} function instead.

If \code{plot} is set to \code{TRUE} and the method is dual pollutant models, then the user will have to choose which pollutant should be plotted. It is not possible to plot all of them at once. A handy menu is provided for selection.

A 2-level interaction term can be supplied in \code{interaction}. The interaction effect will be estimated as well as the marginal effects. The argument \code{interaction} must be a factor or it will be coerced to a factor. Interaction estimation is available for single lag method only.

Over-dispersed models can be fitted by \emph{quasi-likelihood} if \code{overdispersion=TRUE}. One should get larger confidence intervals when this options is set under over-dispersed data.
}
\value{
The function invisibly returns an array of matrices with the exposures on the rows ant the statistics on the columns.
}
\references{
Schwartz, J., Spix, C., Touloumi, G. et al. (1996) Methodological issues in studies of air pollution and daily counts of deaths or hospital admissions. \emph{J Epidemiol. Community Health} \bold{50 (suppl 1)}, S12--S18.

Schwartz, J. (2000) The distributed lag between air pollution and daily deaths. \emph{Epidemiology} \bold{11(3)}, 320--326.

McGullagh, P., Nelder, J. A. (1989) \emph{Generalized linear models}. Chapman and Hall.

Hastie, T., Tibshirani, R. (1990) \emph{Generalized additive models}. 2 ed. Chapman and Hall.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

}

\seealso{\code{\link{fit_core}},\code{\link{print_risk}},\code{\link{plot_risk}}}
\examples{
data(admrio)
setup(admrio,"date")
f <- resp5~s(time)+weekdays+s(tmpmax)+s(humid)
m <- fit_core(f)
## single lag effect estimation
r1 <- estimate_risks(m,c("pm10","so2"),digits=3,labels=c("PM10","SO2"),method="singlelag",
	lag.struc=list(l=0:2,ma=1:5))

## pdlm effect estimation
r2 <- estimate_risks(m,c("pm10","so2"),digits=3,labels=c("PM10","SO2"),method="pdlm",
	pdlm.struc=list(l=5,deg=c(2,2)))

## dual pollutant model (it is commented in order to not run during check)
## r3 <- estimate_risks(m,c("pm10","so2","co"),digits=3,labels=c("PM10","SO2","CO"),
##	method="dual",lag.struc=list(l=0:2,ma=1:5))
}

\keyword{models}
\keyword{regression}
